function out = saSPT_bootstrap_CI(tracksByFOV, varargin)
% SASPT_BOOTSTRAP_CI  Cluster bootstrap CIs for the saSPT posterior.
%
%   out = saSPT_bootstrap_CI(tracksByFOV, 'Name', value, ...)
%
%   Resamples FOVs (wrappers) with replacement, recomputes the saSPT
%   posterior, and returns per-D-bin percentile confidence intervals on the
%   marginalized occupancy `sum(post_occ, 2)`.  Captures between-cell
%   heterogeneity, which is the dominant source of variance for SMT data.
%
%   The likelihood `log_L` is computed once on every subtrack across all
%   FOVs and cached; each bootstrap iteration re-slices it and re-runs only
%   the (cheap) variational posterior update.  This makes the bootstrap
%   roughly free relative to the up-front likelihood evaluation.
%
%   Inputs
%     tracksByFOV : one of the following
%                   - TrajectoryCollection  (uses culled tracks per wrapper)
%                   - cell of cells: {fov1_tracks, fov2_tracks, ...}
%                     where each fovK_tracks is itself a cell array of
%                     [x, y, frame, ...] matrices (µm, frame).
%
%   Name-value parameters (saSPT model)
%     'Condition'      []                              filter (TC input only)
%     'SplitSize'      7
%     'D'              logspace(log10(1e-3), log10(20), 64)   (µm² s⁻¹)
%     'LocError'       linspace(0.01, 0.1, 24)                (µm)
%     'FrameInterval'  []  (required if cell input; auto for TC)
%     'ExposureTime'   []  defaults to FrameInterval; 0 disables blur
%     'ShutterR'       1/6
%
%   Name-value parameters (bootstrap)
%     'NumBootstraps'  200
%     'AlphaCI'        0.05    (95% percentile interval)
%     'Mode'           'cluster'   'cluster' = resample with replacement;
%                                  'jackknife' = leave-one-cluster-out
%     'ClusterBy'      'fov'       resampling unit; one of:
%                                    'fov'   = one cluster per FOV (default).
%                                    'cell'  = one cluster per (FOV, ROI)
%                                              pair.  Requires TC input
%                                              and culled tracks.  Use when
%                                              cells are similarly sized.
%                                    'chunk' = random chunks of ChunkSize
%                                              subtracks each.  Use when
%                                              both FOVs and cells are
%                                              imbalanced (e.g. cells with
%                                              1 vs 1000 tracks).  Equal-
%                                              sized chunks make the
%                                              bootstrap mean concentrate
%                                              tightly around the point
%                                              estimate, at the cost of
%                                              ignoring biological
%                                              clustering — interpret as a
%                                              variance estimate of the
%                                              variational posterior, not
%                                              of cell-to-cell biology.
%     'ChunkSize'      5000        subtracks per chunk (chunk mode only)
%     'MinSubtracksPerCluster'     drop clusters with fewer than this many
%                                  subtracks before bootstrapping.  Default:
%                                  50 for ClusterBy='cell', 0 otherwise.
%                                  Ignored in chunk mode (chunks are
%                                  uniformly sized).
%     'Seed'           []
%     'Verbose'        true
%     'Parallel'       false   use parfor for the bootstrap loop; auto-launches
%                              a parpool if none is open
%
%   Output struct fields
%     D, LocError, FrameInterval, ExposureTime, ShutterR, SplitSize, Mode
%     NumFOVs, NumSubtracksPerFOV
%     post_marg_D_full   LD x 1   posterior on the full (un-resampled) data
%     post_marg_D_boot   LD x B   per-iteration marginal posteriors
%     post_marg_D_mean   LD x 1   mean across iterations
%     post_marg_D_std    LD x 1   std (or jackknife SE)
%     post_marg_D_CI     LD x 2   [lower, upper] percentile CIs

    p = inputParser;
    addParameter(p, 'Condition',     [], @(x) ischar(x)||isstring(x)||iscell(x));
    addParameter(p, 'SplitSize',      7, @(x) isnumeric(x)&&isscalar(x)&&x>=3);
    addParameter(p, 'D',             logspace(log10(1e-3), log10(20), 64), @isnumeric);
    addParameter(p, 'LocError',      linspace(0.01, 0.1, 24),               @isnumeric);
    addParameter(p, 'FrameInterval', [], @(x) isempty(x)||(isscalar(x)&&x>0));
    addParameter(p, 'ExposureTime',  [], @(x) isempty(x)||(isscalar(x)&&x>=0));
    addParameter(p, 'ShutterR',     1/6, @(x) isscalar(x)&&x>=0);
    addParameter(p, 'NumBootstraps', 200, @(x) isscalar(x)&&x>=2);
    addParameter(p, 'AlphaCI',      0.05, @(x) isscalar(x)&&x>0&&x<1);
    addParameter(p, 'Mode',     'cluster', @(s) any(strcmpi(s,{'cluster','jackknife'})));
    addParameter(p, 'ClusterBy',   'fov', @(s) any(strcmpi(s,{'fov','cell','chunk'})));
    addParameter(p, 'ChunkSize',   5000,  @(x) isnumeric(x)&&isscalar(x)&&x>=100);
    addParameter(p, 'MinSubtracksPerCluster', [], @(x) isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=0));
    addParameter(p, 'Seed',          [], @(x) isempty(x)||isnumeric(x));
    addParameter(p, 'Verbose',     true, @islogical);
    addParameter(p, 'Parallel',   false, @islogical);
    parse(p, varargin{:});
    o = p.Results;

    %% Resolve input -> {fov1_tracks, ...}
    if isa(tracksByFOV, 'TrajectoryCollection')
        tc = tracksByFOV;
        if ~isempty(o.Condition)
            mask = strcmp(tc.Metadata.Condition, char(o.Condition));
            wrappers = tc.Wrappers(mask);
        else
            wrappers = tc.Wrappers;
        end
        fovTracks = cellfun(@(w) w.getCulledTracks(), wrappers, 'UniformOutput', false);
        if isempty(o.FrameInterval)
            o.FrameInterval = tc.getFrameIntervalForCondition(o.Condition);
        end
    elseif iscell(tracksByFOV) && ~isempty(tracksByFOV) && iscell(tracksByFOV{1})
        fovTracks = tracksByFOV(:)';
    else
        error('tracksByFOV must be a TrajectoryCollection or a cell of cell arrays.');
    end

    if isempty(o.FrameInterval) || isnan(o.FrameInterval)
        error('FrameInterval must be supplied (or set on wrappers).');
    end
    if isempty(o.ExposureTime)
        o.ExposureTime = o.FrameInterval;
    end

    nFOV = numel(fovTracks);
    if nFOV < 2
        error('Need at least 2 FOVs; got %d.', nFOV);
    end

    %% Build cluster track-groups
    %  ClusterBy='fov'  : one cluster per FOV
    %  ClusterBy='cell' : one cluster per (FOV, ROI) pair.  ROI IDs are
    %                     unique within a FOV but repeat across FOVs, so the
    %                     global cluster key is (FOV index, ROI ID).
    if strcmpi(o.ClusterBy, 'cell')
        if ~isa(tracksByFOV, 'TrajectoryCollection')
            error(['ClusterBy=''cell'' requires a TrajectoryCollection input ' ...
                   '(ROI IDs are not available from raw cell-of-cell input).']);
        end
        clusterTracks = {};
        clusterLabel  = {};
        for k = 1:nFOV
            wTracks = fovTracks{k};
            if isempty(wTracks), continue; end
            badShape = cellfun(@(t) size(t,2) < 5, wTracks);
            if any(badShape)
                error(['ClusterBy=''cell'' requires culled tracks (Nx5) with ' ...
                       'ROI ID in column 5; FOV %d has tracks without it.'], k);
            end
            roiVec = cellfun(@(t) t(1,5), wTracks);
            uROI   = unique(roiVec);
            for r = 1:numel(uROI)
                mask = roiVec == uROI(r);
                clusterTracks{end+1} = wTracks(mask);                   %#ok<AGROW>
                clusterLabel{end+1}  = sprintf('FOV%d/ROI%d', k, uROI(r)); %#ok<AGROW>
            end
        end
    elseif strcmpi(o.ClusterBy, 'chunk')
        % For chunk mode, pool everything into a single "cluster" so the
        % likelihood is computed once on the full subtrack pool.  Random
        % chunks of size ChunkSize are formed AFTER the likelihood, by
        % re-partitioning the (D, sigma, subtrack) array along the third
        % dimension.  See the post-likelihood block below.
        % Some wrappers return culled tracks as a column cell, others as a
        % row cell; horzcat requires uniform shape, so reshape to row first.
        fovRow = cellfun(@(c) reshape(c, 1, []), fovTracks, ...
                         'UniformOutput', false);
        clusterTracks = { horzcat(fovRow{:}) };
        clusterLabel  = {'all_pooled'};
    else
        clusterTracks = fovTracks;
        clusterLabel  = arrayfun(@(k) sprintf('FOV%d', k), 1:nFOV, 'UniformOutput', false);
    end

    nClust = numel(clusterTracks);
    % For chunk mode, the single pooled "cluster" is an intermediate; the
    % real cluster count is set after the post-likelihood re-partition below.
    if nClust < 2 && ~strcmpi(o.ClusterBy, 'chunk')
        error('Need at least 2 clusters for bootstrap; got %d (ClusterBy=''%s'').', ...
              nClust, lower(o.ClusterBy));
    end

    %% Split + likelihood once per cluster
    if o.Verbose
        fprintf('saSPT bootstrap: %d clusters (ClusterBy=%s), splitSize=%d, |D|=%d, |sigma|=%d\n', ...
            nClust, lower(o.ClusterBy), o.SplitSize, numel(o.D), numel(o.LocError));
    end

    fovLogL    = cell(1, nClust);
    fovJumps   = cell(1, nClust);
    nSubPerFOV = zeros(1, nClust);

    if o.Verbose
        if strcmpi(o.ClusterBy, 'chunk')
            fprintf('  computing likelihood on the full pooled subtrack set...\n');
        else
            fprintf('  computing per-cluster likelihoods (%d clusters)...\n', nClust);
        end
        t_lik = tic;
        likStep = max(1, round(nClust/20));
    end
    for k = 1:nClust
        sub = saSPT_splitTracks(clusterTracks{k}, o.SplitSize);
        nSubPerFOV(k) = numel(sub);
        if isempty(sub)
            fovLogL{k}  = zeros(numel(o.D), numel(o.LocError), 0);
            fovJumps{k} = zeros(0, 1);
            continue;
        end
        [fovLogL{k}, fovJumps{k}] = saSPT_calc_likelihood( ...
            sub, o.D, o.LocError, o.FrameInterval, o.ExposureTime, o.ShutterR);
        if o.Verbose
            if nClust <= 50
                fprintf('    %s: %d subtracks\n', clusterLabel{k}, nSubPerFOV(k));
            elseif mod(k, likStep) == 0 || k == nClust
                elapsed = toc(t_lik);
                eta     = elapsed * (nClust - k) / max(k, 1);
                fprintf('    %4d/%d clusters | elapsed=%5.1fs  ETA=%5.1fs\n', ...
                    k, nClust, elapsed, eta);
            end
        end
    end
    if o.Verbose
        fprintf('  per-cluster likelihoods done (%.1fs).\n', toc(t_lik));
    end

    %% Drop clusters that are too small to contribute meaningfully.
    %  Default threshold: 50 subtracks for ClusterBy='cell' (where many
    %  cells contain only a handful of tracks), 0 otherwise.  User can
    %  override either way via 'MinSubtracksPerCluster'.
    minSubs = o.MinSubtracksPerCluster;
    if isempty(minSubs)
        if strcmpi(o.ClusterBy, 'cell'), minSubs = 50; else, minSubs = 0; end
    end
    if minSubs > 0 && ~strcmpi(o.ClusterBy, 'chunk')
        keep = nSubPerFOV >= minSubs;
        nDropped = sum(~keep);
        if sum(keep) < 2
            error(['MinSubtracksPerCluster=%d leaves only %d clusters; ' ...
                   'lower the threshold.'], minSubs, sum(keep));
        end
        if o.Verbose && nDropped > 0
            droppedSubs = sum(nSubPerFOV(~keep));
            fprintf(['  dropping %d/%d clusters with <%d subtracks ' ...
                     '(%d subtracks total, %.1f%% of pool).\n'], ...
                nDropped, nClust, minSubs, droppedSubs, ...
                100*droppedSubs/max(1,sum(nSubPerFOV)));
        end
        fovLogL      = fovLogL(keep);
        fovJumps     = fovJumps(keep);
        nSubPerFOV   = nSubPerFOV(keep);
        clusterLabel = clusterLabel(keep);
        nClust       = numel(fovLogL);
    end

    %% Chunk mode: re-partition the pooled likelihood into random chunks
    %  of size ChunkSize.  Each subtrack is independently assigned to a
    %  chunk, so chunk size is uniform (last chunk may be smaller).
    if strcmpi(o.ClusterBy, 'chunk')
        L_pool = fovLogL{1};       % [|D|, |sigma|, N]
        j_pool = fovJumps{1};      % [N, 1]
        N      = size(L_pool, 3);
        if N < 2*o.ChunkSize
            error(['ClusterBy=''chunk'' with ChunkSize=%d needs at least %d ' ...
                   'subtracks for >=2 chunks; got %d.  Lower ChunkSize.'], ...
                   o.ChunkSize, 2*o.ChunkSize, N);
        end
        if ~isempty(o.Seed), rng(o.Seed); end
        perm     = randperm(N);
        nChunks  = ceil(N / o.ChunkSize);
        fovLogL    = cell(1, nChunks);
        fovJumps   = cell(1, nChunks);
        nSubPerFOV = zeros(1, nChunks);
        clusterLabel = cell(1, nChunks);
        for c = 1:nChunks
            lo = (c-1)*o.ChunkSize + 1;
            hi = min(c*o.ChunkSize, N);
            idx = perm(lo:hi);
            fovLogL{c}     = L_pool(:, :, idx);
            fovJumps{c}    = j_pool(idx);
            nSubPerFOV(c)  = numel(idx);
            clusterLabel{c} = sprintf('Chunk%d', c);
        end
        nClust = nChunks;
        if o.Verbose
            fprintf('  partitioned %d subtracks into %d chunks of ~%d (last=%d).\n', ...
                N, nChunks, o.ChunkSize, nSubPerFOV(end));
        end
    end
    nNonEmpty = sum(nSubPerFOV > 0);
    nSubTotal = sum(nSubPerFOV);
    if o.Verbose
        if nClust > 50
            fprintf('  per-cluster subtrack counts: median=%d, min=%d, max=%d (max/min=%.1fx)\n', ...
                round(median(nSubPerFOV(nSubPerFOV>0))), ...
                min(nSubPerFOV(nSubPerFOV>0)), max(nSubPerFOV), ...
                max(nSubPerFOV) / max(1, min(nSubPerFOV(nSubPerFOV>0))));
        end
        fprintf('  total: %d subtracks across %d non-empty clusters (of %d).\n', ...
            nSubTotal, nNonEmpty, nClust);
    end

    if all(nSubPerFOV == 0)
        error('No subtracks produced from any cluster.');
    end

    %% Posterior on the full pool (point estimate)
    if o.Verbose
        fprintf('  computing point-estimate posterior on full pool...\n');
        t_full = tic;
    end
    log_L_all = cat(3, fovLogL{:});
    jumps_all = vertcat(fovJumps{:});
    [~, post_full] = saSPT_compute_posterior(log_L_all, jumps_all, o.D, o.LocError);
    post_marg_full = sum(post_full, 2);
    if o.Verbose
        fprintf('  point estimate done (%.1fs).\n', toc(t_full));
    end

    %% Bootstrap iterations
    LD = numel(o.D);
    switch lower(o.Mode)
        case 'cluster'
            B = o.NumBootstraps;
            if ~isempty(o.Seed), rng(o.Seed); end
            idxMatrix = randi(nClust, B, nClust);
        case 'jackknife'
            B = nClust;
            idxMatrix = zeros(B, nClust - 1);
            for i = 1:B
                idxMatrix(i, :) = setdiff(1:nClust, i);
            end
    end

    post_marg_boot = zeros(LD, B);

    %% Launch parallel pool if requested
    useParallel = o.Parallel;
    if useParallel
        if isempty(ver('parallel'))
            warning('Parallel Computing Toolbox not available; falling back to serial.');
            useParallel = false;
        elseif isempty(gcp('nocreate'))
            parpool;
        end
    end

    if o.Verbose
        modeStr = sprintf('%s, ClusterBy=%s, %d clusters', ...
            lower(o.Mode), lower(o.ClusterBy), nClust);
        if useParallel
            fprintf('Running %d bootstrap iterations (%s, parallel)...\n', B, modeStr);
        else
            fprintf('Running %d bootstrap iterations (%s, serial)...\n', B, modeStr);
        end
    end

    % Broadcast vars for parfor
    fovLogL_b   = fovLogL;
    fovJumps_b  = fovJumps;
    Dgrid       = o.D;
    sigGrid     = o.LocError;
    verbose_b   = o.Verbose && ~useParallel;
    progressStep = max(1, round(B/20));   % print ~20 updates total
    boot_nsubs   = zeros(1, B);

    t_boot = tic;
    if useParallel
        % DataQueue for progress
        dq = parallel.pool.DataQueue;
        afterEach(dq, @(payload) print_progress_par(payload, B, progressStep, t_boot));
        parfor b = 1:B
            sel = idxMatrix(b, :);
            log_L_b = cat(3, fovLogL_b{sel});
            jumps_b = vertcat(fovJumps_b{sel});
            ns = numel(jumps_b);
            boot_nsubs(b) = ns;
            if isempty(jumps_b)
                post_marg_boot(:, b) = NaN;
            else
                post_b = quiet_posterior(log_L_b, jumps_b, Dgrid, sigGrid);
                post_marg_boot(:, b) = sum(post_b, 2);
            end
            send(dq, [b, ns]);
        end
    else
        for b = 1:B
            sel = idxMatrix(b, :);
            log_L_b = cat(3, fovLogL_b{sel});
            jumps_b = vertcat(fovJumps_b{sel});
            ns = numel(jumps_b);
            boot_nsubs(b) = ns;
            if isempty(jumps_b)
                post_marg_boot(:, b) = NaN;
                continue;
            end
            post_b = quiet_posterior(log_L_b, jumps_b, Dgrid, sigGrid);
            post_marg_boot(:, b) = sum(post_b, 2);
            if verbose_b && (mod(b, progressStep) == 0 || b == B)
                elapsed = toc(t_boot);
                eta     = elapsed * (B - b) / max(b, 1);
                fprintf('  iter %4d/%d  | nSub=%d  | elapsed=%5.1fs  ETA=%5.1fs\n', ...
                    b, B, ns, elapsed, eta);
            end
        end
    end
    if o.Verbose
        fprintf('Bootstrap done (%.1fs total). Resample size: median=%d, range=[%d, %d].\n', ...
            toc(t_boot), round(median(boot_nsubs)), min(boot_nsubs), max(boot_nsubs));
    end

    %% Aggregate
    post_marg_mean = mean(post_marg_boot, 2, 'omitnan');
    switch lower(o.Mode)
        case 'cluster'
            post_marg_std = std(post_marg_boot, 0, 2, 'omitnan');
            qLo = o.AlphaCI/2;  qHi = 1 - o.AlphaCI/2;
            post_marg_CI = quantile(post_marg_boot, [qLo qHi], 2);
        case 'jackknife'
            % Jackknife SE = sqrt((N-1)/N * sum((xi - mean)^2))
            d = post_marg_boot - post_marg_mean;
            post_marg_std = sqrt((nClust-1)/nClust * sum(d.^2, 2, 'omitnan'));
            % Approx normal CI from jackknife SE
            z = norminv(1 - o.AlphaCI/2);
            post_marg_CI = [post_marg_full - z*post_marg_std, post_marg_full + z*post_marg_std];
    end

    out = struct( ...
        'D',                      o.D, ...
        'LocError',               o.LocError, ...
        'FrameInterval',          o.FrameInterval, ...
        'ExposureTime',           o.ExposureTime, ...
        'ShutterR',               o.ShutterR, ...
        'SplitSize',              o.SplitSize, ...
        'Condition',              o.Condition, ...
        'Mode',                   lower(o.Mode), ...
        'ClusterBy',              lower(o.ClusterBy), ...
        'NumFOVs',                nFOV, ...
        'NumClusters',            nClust, ...
        'ClusterLabels',         {clusterLabel}, ...
        'NumSubtracksPerCluster', nSubPerFOV, ...
        'NumSubtracksPerFOV',     nSubPerFOV, ...   % alias kept for back-compat
        'opts',                   o, ...
        'post_marg_D_full',       post_marg_full, ...
        'post_marg_D_boot',       post_marg_boot, ...
        'post_marg_D_mean',       post_marg_mean, ...
        'post_marg_D_std',        post_marg_std, ...
        'post_marg_D_CI',         post_marg_CI);
end

function post_b = quiet_posterior(log_L_b, jumps_b, Dgrid, sigGrid)
    % Wrap evalc inside a function so its workspace mutation does not
    % violate parfor transparency rules.
    evalc('[~, post_b] = saSPT_compute_posterior(log_L_b, jumps_b, Dgrid, sigGrid);');
end

function print_progress_par(payload, B, step, t_boot)
    % payload = [iteration_index, n_subtracks_in_resample]
    persistent count last_ns
    if isempty(count),   count = 0;   end
    if isempty(last_ns), last_ns = 0; end
    count   = count + 1;
    last_ns = payload(2);
    if mod(count, step) == 0 || count == B
        elapsed = toc(t_boot);
        eta     = elapsed * (B - count) / max(count, 1);
        fprintf('  iter %4d/%d  | nSub=%d  | elapsed=%5.1fs  ETA=%5.1fs\n', ...
            count, B, last_ns, elapsed, eta);
    end
    if count == B, count = 0; last_ns = 0; end
end
