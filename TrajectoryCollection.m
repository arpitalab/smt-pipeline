classdef TrajectoryCollection < handle
    % TRAJECTORYCOLLECTION - Aggregates multiple TrajectoryWrapper instances for
    %                        multi-file / multi-condition SMT analysis.
    %
    % Supports two input paths:
    %   CSV path  : tc.addFromFile('locs.csv', ...)
    %   SMD path  : tc.addFromSMD(smd, ...)
    %
    % Both paths produce identical TrajectoryWrapper objects that flow into
    % the same MSD / Bayesian / lifetime analysis methods.
    %
    % Example (mixed sources):
    %   tc = TrajectoryCollection();
    %   tc.ReadParams('experiment.toml');
    %   tc.addFromFile('s001.csv', 'FileID','1', 'Condition','Control');
    %   tc.addFromSMD(smd,         'FileID','2', 'Condition','Treated');
    %   tc.getBayesianDiffusivity('Condition','Control');
    %   tc.summary();

    properties (SetAccess = private)
        Wrappers           = {}           % Cell array of TrajectoryWrapper instances
        Metadata                          % Table: FileID, FilePath, Condition
        AllRawTracks       = {}           % Lazy: aggregated raw tracks
        AllCulledTracks    = {}           % Lazy: aggregated culled tracks
        AllRelativeTracks  = {}           % Lazy: aggregated relative tracks
        RLResults          = struct()
        MSDResults         = struct()
        pEMResults         = struct()
        Params             = struct()
        pEMParams          = struct()
        LifetimeResults    = struct()
    end

    properties (Access = private)
        IsAllRawCollected      = false
        IsAllCulledCollected   = false
        IsAllRelativeCollected = false
        IsRLComputed           = false
        IsMSDComputed          = false
        IspEMComputed          = false
        IsBootstrapComputed    = false
        IsLifetimeComputed     = false
        pEMBootstrapInputs     = struct()  % deltaX, rawResults, trackInfo, params
    end

    methods
        %% Constructor
        function obj = TrajectoryCollection(varargin)
            p = inputParser;
            addParameter(p, 'Condition',  '');
            addParameter(p, 'Parameters', struct());
            parse(p, varargin{:});

            obj.Metadata = table('Size', [0 3], ...
                'VariableTypes', {'string','string','string'}, ...
                'VariableNames', {'FileID','FilePath','Condition'});

            if ~isempty(p.Results.Condition)
                obj.Metadata.Condition(:) = p.Results.Condition;
            end
            if ~isempty(fieldnames(p.Results.Parameters))
                obj.addGlobalParameters(p.Results.Parameters);
            end
        end

        %% Add from CSV file
        function obj = addFromFile(obj, csvFile, varargin)
            % ADDFROMFILE  Create a TrajectoryWrapper from a CSV file, cull, and add.
            %
            %   tc.addFromFile(csvFile)
            %   tc.addFromFile(csvFile, 'Condition', 'Control', ...
            %                           'PixelSize', 0.165, ...
            %                           'FrameInterval', 0.02, ...
            %                           'Parameters', params)
            %   tc.addFromFile(csvFile, 'ParamsFile', 'experiment.toml', ...)
            %
            %   PixelSize and FrameInterval can be supplied directly or read
            %   from a TOML file via 'ParamsFile'.  Direct values take
            %   precedence over TOML values.  PixelSize must be known before
            %   readData() runs (it scales x,y from pixels to µm on load).

            p = inputParser;
            addParameter(p, 'FileID',        '');
            addParameter(p, 'Condition',      '');
            addParameter(p, 'Parameters',     struct());
            addParameter(p, 'ParamsFile',     '');       % path to TOML file
            addParameter(p, 'PixelSize',      NaN);      % µm/px; NaN = use TOML or TW default
            addParameter(p, 'FrameInterval',  NaN);      % s;     NaN = use TOML or leave unset
            parse(p, varargin{:});

            fileID = p.Results.FileID;
            if isempty(fileID)
                [~, fileID, ~] = fileparts(csvFile);
            end

            tw = TrajectoryWrapper('Parameters', p.Results.Parameters);

            % Load TOML first (lowest precedence)
            if ~isempty(p.Results.ParamsFile)
                tw.readParams(p.Results.ParamsFile);
            end

            % Direct arguments override TOML values
            if ~isnan(p.Results.PixelSize)
                tw.setPixelSize(p.Results.PixelSize);
            end
            if ~isnan(p.Results.FrameInterval)
                tw.setFrameInterval(p.Results.FrameInterval);
            end

            tw.readData(csvFile);
            try
                tw.loadNucleusMask();
                tw.cull_tracks();
                tw.remove_relative_motion();
            catch
                warning('no nucleus file or cullable tracks');
            end

            obj = obj.addWrapper(tw, 'FilePath', csvFile, ...
                'FileID', fileID, 'Condition', p.Results.Condition);
        end

        %% Add from a processed SMD object
        function obj = addFromSMD(obj, smd, varargin)
            % ADDFROMSMD  Create a TrajectoryWrapper from a processed SMD object and add.
            %
            %   tc.addFromSMD(smd)
            %   tc.addFromSMD(smd, 'FileID','run1', 'Condition','Control')
            %
            %   The SMD object must already be processed (localize → track →
            %   get_roi → cull_tracks). No pipeline steps are run here.

            p = inputParser;
            addParameter(p, 'FileID',    '');
            addParameter(p, 'Condition', '');
            parse(p, varargin{:});

            fileID = p.Results.FileID;
            if isempty(fileID)
                [~, fileID, ~] = fileparts(smd.fname);
            end

            tw = TrajectoryWrapper.fromSMD(smd);

            obj = obj.addWrapper(tw, 'FilePath', smd.fname, ...
                'FileID', fileID, 'Condition', p.Results.Condition);
        end

        %% Add existing TrajectoryWrapper
        function obj = addWrapper(obj, tw, varargin)
            if ~isa(tw, 'TrajectoryWrapper')
                error('Input must be a TrajectoryWrapper instance.');
            end
            p = inputParser;
            addParameter(p, 'FileID',   num2str(length(obj.Wrappers) + 1));
            addParameter(p, 'FilePath', '');
            addParameter(p, 'Condition','');
            parse(p, varargin{:});

            obj.Wrappers{end+1} = tw;
            obj.Metadata = [obj.Metadata; {p.Results.FileID, p.Results.FilePath, p.Results.Condition}];
            obj.invalidateCollections();
        end

        %% Batch load from file lists
        function obj = loadFromFiles(obj, csvFiles, maskFiles, checkFunc, conditions)
            if length(csvFiles) ~= length(maskFiles)
                error('CSV and mask file lists must match in length.');
            end
            if nargin < 5 || isempty(conditions)
                conditions = repmat({''}, length(csvFiles), 1);
            elseif ischar(conditions) || isstring(conditions)
                conditions = repmat({conditions}, length(csvFiles), 1);
            end
            for i = 1:length(csvFiles)
                obj.addFromFile(csvFiles{i}, maskFiles{i}, checkFunc, ...
                    'FileID', num2str(i), 'Condition', conditions{i});
            end
        end

        %% Read experiment params from TOML
        function obj = ReadParams(obj, fileID)
            % READPARAMS  Load TOML parameters into Params and propagate to wrappers.
            %
            %   Recognized top-level TOML keys (propagated to wrappers):
            %     pixel_size      → sets PixelSize on CSV-sourced wrappers
            %     frame_interval  → sets FrameInterval on wrappers where it is NaN
            %   All keys are also stored in obj.Params for analysis code to read.

            config     = toml.read(fileID);
            pars       = toml.map_to_struct(config);
            obj.Params = pars;

            % Propagate pixel_size / frame_interval to existing wrappers
            for w = 1:length(obj.Wrappers)
                tw = obj.Wrappers{w};
                if isfield(pars, 'pixel_size') && strcmp(tw.getSourceType(), 'csv')
                    % Only override for CSV wrappers; SMD wrappers derive px from camera cal
                    tw.setPixelSize(pars.pixel_size);
                end
                if isfield(pars, 'frame_interval') && isnan(tw.getFrameInterval())
                    tw.setFrameInterval(pars.frame_interval);
                end
            end
            obj.invalidateCollections();
        end

        %% Set frame interval on all wrappers
        function obj = setFrameIntervalAll(obj, dt)
            % SETFRAMEINTERVALALL  Set FrameInterval (s) on every wrapper in the collection.
            %
            %   tc.setFrameIntervalAll(0.02)
            %
            %   Overrides existing values on all wrappers (including SMD-sourced ones).
            %   Call this after all addFromFile / addFromSMD calls are complete.
            for w = 1:length(obj.Wrappers)
                obj.Wrappers{w}.setFrameInterval(dt);
            end
            obj.invalidateCollections();
        end

        %% Lazy accessors
        function tracks = getAllRawTracks(obj)
            if ~obj.IsAllRawCollected, obj.collectAllTracks(); end
            tracks = obj.AllRawTracks;
        end

        function tracks = getAllCulledTracks(obj)
            if ~obj.IsAllCulledCollected, obj.collectAllTracks(); end
            tracks = obj.AllCulledTracks;
        end

        function tracks = getAllRelativeTracks(obj)
            if ~obj.IsAllRelativeCollected, obj.collectAllTracks(); end
            tracks = obj.AllRelativeTracks;
        end

        %% RL decomposition (Richardson-Lucy MSD distribution)
        function results = getRLDecomposition(obj, varargin)
            % GETRLDECOMPOSITION  Richardson-Lucy decomposition of the MSD distribution.
            %
            %   results = tc.getRLDecomposition('LagTime', 4, 'String', 'H2B')
            %
            %   Runs RL_HILO on culled tracks (or a condition subset) and stores
            %   results in tc.RLResults.  Re-call with the same arguments returns
            %   cached results; use 'ForceRecompute' to re-run.
            %
            %   Options:
            %     'LagTime'        - frame lag for van Hove computation (default 4)
            %     'String'         - label for figure titles (default 'Not specified')
            %     'Condition'      - restrict to one condition (default: all)
            %     'ForceRecompute' - logical; default false

            if ~obj.IsRLComputed
                p = inputParser;
                addParameter(p, 'LagTime',        4);
                addParameter(p, 'String',          'Not specified');
                addParameter(p, 'Condition',       [], @(x) ischar(x)||isstring(x)||iscell(x));
                addParameter(p, 'ForceRecompute',  false, @islogical);
                parse(p, varargin{:});

                if isempty(p.Results.Condition)
                    tracks = obj.getAllCulledTracks();
                else
                    tracks = obj.getTracksByCondition(p.Results.Condition);
                end

                if isempty(tracks)
                    error('No tracks available for RL decomposition.');
                end

                X = cell(1, numel(tracks));
                for itrack = 1:numel(tracks)
                    X{itrack} = tracks{itrack}(:, 1:2);   % already in µm
                end
                tracklen = cellfun('size', tracks, 1);
                idx      = find(tracklen > 6);
                tmp{1}   = X(idx);
                [M, ~, computed_quantities] = RL_HILO(tmp, p.Results.String, p.Results.LagTime);
                obj.RLResults.M                   = M;
                obj.RLResults.computed_quantities = computed_quantities;
                obj.IsRLComputed = true;
            end
            results = obj.RLResults;
        end

        %% Ensemble MSD with bootstrap fitting
        function results = getMSD(obj, varargin)
            % GETMSD  Ensemble-averaged MSD with bootstrap power-law fits.
            %
            %   results = tc.getMSD()
            %   results = tc.getMSD('TrackType',    'culled')   % 'culled'|'relative'|'raw'
            %   results = tc.getMSD('MaxLag',       25)         % max lag in frames
            %   results = tc.getMSD('NumBootstrap', 50)         % bootstrap samples
            %   results = tc.getMSD('ExposureTime', 0.01)       % camera exposure (s)
            %   results = tc.getMSD('Condition',    'Control')
            %   results = tc.getMSD('ForceRecompute', true)
            %
            %   ExposureTime defaults to Params.exposure_time_s if set via TOML,
            %   otherwise to the frame interval (frac = 1, i.e. stroboscopic).
            %
            %   Results struct fields:
            %     .ensemble_mean  - [MaxLag x 1] ensemble-averaged MSD (µm²)
            %     .lag_axis       - [MaxLag x 1] lag times in seconds
            %     .lag_frames     - [MaxLag x 1] lag indices (1 … MaxLag)
            %     .boot_means     - [MaxLag x NumBootstrap] per-bootstrap means
            %     .boot_fits      - [NumBootstrap x 3] fit params [G, sigma_sq, alpha]
            %     .dt             - frame interval (s)
            %     .frac           - exposure / frame_interval ratio
            %     .track_type     - which track set was used

            if ~obj.IsMSDComputed
                p = inputParser;
                addParameter(p, 'TrackType',      'culled', @(x) ischar(x)||isstring(x));
                addParameter(p, 'MaxLag',         25,       @isnumeric);
                addParameter(p, 'NumBootstrap',   50,       @isnumeric);
                addParameter(p, 'ExposureTime',   NaN,      @isnumeric);
                addParameter(p, 'Condition',      [],       @(x) ischar(x)||isstring(x)||iscell(x));
                addParameter(p, 'ForceRecompute', false,    @islogical);
                parse(p, varargin{:});
                o = p.Results;

                % --- Get tracks -------------------------------------------
                trackType = lower(char(o.TrackType));
                switch trackType
                    case 'culled'
                        if isempty(o.Condition)
                            tracks = obj.getAllCulledTracks();
                        else
                            tracks = obj.getTracksByCondition(o.Condition);
                        end
                    case 'relative'
                        if isempty(o.Condition)
                            tracks = obj.getAllRelativeTracks();
                        else
                            tracks = obj.getRelativeTracksByCondition(o.Condition);
                        end
                    case 'raw'
                        if isempty(o.Condition)
                            tracks = obj.getAllRawTracks();
                        else
                            tracks = obj.getTracksByCondition(o.Condition, 'UseRaw', true);
                        end
                    otherwise
                        error('getMSD: TrackType must be ''culled'', ''relative'', or ''raw''.');
                end

                if isempty(tracks)
                    error('getMSD: no tracks available (TrackType=''%s'').', trackType);
                end

                % --- Physical parameters ----------------------------------
                dt = obj.getFrameIntervalForCondition(o.Condition);

                expT = o.ExposureTime;
                if isnan(expT) && isfield(obj.Params, 'exposure_time_s')
                    expT = obj.Params.exposure_time_s;
                end
                if isnan(expT)
                    frac = 1;   % assume stroboscopic (exposure = frame interval)
                    warning('getMSD:noExposureTime', ...
                        ['ExposureTime not set; assuming frac=1 (stroboscopic). ' ...
                         'Pass ExposureTime or set exposure_time_s in your TOML.']);
                else
                    frac = expT / dt;
                end

                MaxLag = o.MaxLag;

                % --- Gap-fill and per-track MSD ---------------------------
                nTracks      = numel(tracks);
                ensemble_msd  = cell(nTracks, 1);
                ensemble_nlags = cell(nTracks, 1);

                for ii = 1:nTracks
                    t_raw = tracks{ii};
                    x = t_raw(:, 1);
                    y = t_raw(:, 2);
                    if size(t_raw, 2) >= 3
                        t = t_raw(:, 3) - t_raw(1, 3) + 1;   % 1-based frame index
                    else
                        t = (1:length(x))';
                    end
                    [t_u, x_u] = fillGapsWithNaN(t, x);
                    [~,   y_u] = fillGapsWithNaN(t, y);
                    [msd_i, nlags_i] = simple_msd([x_u, y_u], t_u);
                    ensemble_msd{ii}   = msd_i;
                    ensemble_nlags{ii} = nlags_i;
                end

                % --- Ensemble mean ----------------------------------------
                msdlen = cellfun(@length, ensemble_msd);
                N      = max(msdlen);
                e_sum  = zeros(N, 1);   % weighted sum of MSDs
                n_sum  = zeros(N, 1);   % total observation count at each lag

                for ii = 1:nTracks
                    L = msdlen(ii);
                    e_sum(1:L) = e_sum(1:L) + ensemble_msd{ii};
                    n_sum(1:L) = n_sum(1:L) + ensemble_nlags{ii};
                end
                n_sum(n_sum == 0) = NaN;
                ensemble_mean_full = e_sum ./ n_sum;

                % Trim to MaxLag
                nLags         = min(MaxLag, N);
                ensemble_mean = ensemble_mean_full(1:nLags);
                lag_frames    = (1:nLags)';
                lag_axis      = lag_frames * dt;

                % --- Bootstrap -------------------------------------------
                nBoot      = o.NumBootstrap;
                boot_means = zeros(nLags, nBoot);
                boot_fits  = zeros(nBoot, 3);

                lsqOpts = optimoptions('lsqnonlin', ...
                    'MaxFunctionEvaluations', 10000, ...
                    'Display', 'none', ...
                    'Algorithm', 'levenberg-marquardt');

                x0  = [1e-4, 4e-4, 0.5];
                lb  = [1e-5, 1e-6, 0.1];
                ub  = [5,    4e-3, 1.5];

                for isamp = 1:nBoot
                    idx  = datasample(1:nTracks, nTracks);
                    bs_sum  = zeros(N, 1);
                    bs_n    = zeros(N, 1);
                    for jj = 1:nTracks
                        ii = idx(jj);
                        L  = msdlen(ii);
                        bs_sum(1:L) = bs_sum(1:L) + ensemble_msd{ii};
                        bs_n(1:L)   = bs_n(1:L)   + ensemble_nlags{ii};
                    end
                    bs_n(bs_n == 0) = NaN;
                    bm = bs_sum ./ bs_n;
                    boot_means(:, isamp) = bm(1:nLags);

                    bm_fit = bm(1:nLags);
                    valid  = ~isnan(bm_fit) & bm_fit > 0;
                    if sum(valid) >= 3
                        try
                            fp = lsqnonlin( ...
                                @(x) my_fun(x, lag_frames(valid), bm_fit(valid), dt, frac), ...
                                x0, lb, ub, lsqOpts);
                            boot_fits(isamp, :) = fp;
                        catch
                            boot_fits(isamp, :) = NaN;
                        end
                    else
                        boot_fits(isamp, :) = NaN;
                    end
                end

                % --- Store ------------------------------------------------
                obj.MSDResults.ensemble_mean = ensemble_mean;
                obj.MSDResults.lag_axis      = lag_axis;
                obj.MSDResults.lag_frames    = lag_frames;
                obj.MSDResults.boot_means    = boot_means;
                obj.MSDResults.boot_fits     = boot_fits;
                obj.MSDResults.dt            = dt;
                obj.MSDResults.frac          = frac;
                obj.MSDResults.track_type    = trackType;
                obj.MSDResults.n_tracks      = nTracks;
                obj.IsMSDComputed = true;

                fprintf('getMSD: %d tracks, dt=%.4fs, frac=%.3f, MaxLag=%d, %d bootstrap samples.\n', ...
                    nTracks, dt, frac, nLags, nBoot);
            end
            results = obj.MSDResults;
        end

        %% Bayesian diffusivity
        function results = getBayesianDiffusivity(obj, varargin)
            % GETBAYESIANDIFFUSIVITY  Bayesian (pEMv2) diffusivity estimation.
            %
            %   results = tc.getBayesianDiffusivity()
            %   results = tc.getBayesianDiffusivity('Condition','Control')
            %
            % FrameInterval is read per-condition from wrappers.  If any
            % wrapper has FrameInterval=NaN this will error with a clear message.

            p = inputParser;
            addParameter(p, 'Condition',      [], @(x) ischar(x)||isstring(x)||iscell(x));
            addParameter(p, 'ForceRecompute', [], @islogical);
            addParameter(p, 'MinLength',       6, @isnumeric);
            addParameter(p, 'SplitLength',     7, @isnumeric);
            addParameter(p, 'NumFeatures',    [], @isnumeric);  % default: splitLength-1, max 10
            parse(p, varargin{:});
            opts = p.Results;

            isSubset = ~isempty(opts.Condition);

            if isempty(opts.ForceRecompute) && obj.IspEMComputed
                force = false;
            elseif isempty(opts.ForceRecompute) && ~obj.IspEMComputed
                force = true;
            elseif ~isempty(opts.ForceRecompute) && opts.ForceRecompute
                force = true;
            end
            if isSubset
                force = true;
            end

            fprintf('\n=== Bayesian Diffusivity Analysis ===\n');
            fprintf('Timestamp: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

            if isSubset
                condStr = strjoin(cellstr(opts.Condition), ', ');
                fprintf('Restricted to condition(s): %s\n', condStr);
            else
                fprintf('Analyzing full collection (%d FOVs)\n', length(obj.Wrappers));
            end

            if isSubset
                tracks = obj.getTracksByCondition(opts.Condition);
                dt     = obj.getFrameIntervalForCondition(opts.Condition);
            else
                tracks = obj.getAllCulledTracks();
                dt     = obj.getFrameIntervalForCondition([]);
            end

            if isempty(tracks)
                warning('No tracks available for Bayesian analysis.');
                results = struct();
                fprintf('=====================================\n\n');
                return;
            end

            if force
                fprintf('Analyzing (%d tracks)\n', numel(tracks));
                splitLength = opts.SplitLength;
                if isempty(opts.NumFeatures)
                    numFeatures = min(splitLength - 1, 10);
                else
                    numFeatures = opts.NumFeatures;
                end
                dE             = 0.01;
                minStates      = 1;
                maxStates      = 8;
                numReinitialize = 20;
                numPerturb      = 200;
                maxiter         = 10000;
                convergence     = 1e-7;
                lambda          = 0.0;

                k       = 1;
                X       = {};
                trackID = [];
                for itrack = 1:numel(tracks)
                    if size(tracks{itrack}, 1) > opts.MinLength
                        X{k}       = tracks{itrack}(:, 1:2);   % already in µm
                        trackID(k) = k;
                        k = k + 1;
                    end
                end

                [splitX, splitIndex] = SplitTracks(X, splitLength);

                trackInfo.numberOfTracks = length(splitX);
                trackInfo.dimensions     = size(splitX{1}, 2);
                trackInfo.numFeatures    = numFeatures;
                trackInfo.splitLength    = splitLength;
                trackInfo.splitIndex     = splitIndex;
                trackInfo.dt             = dt;
                trackInfo.R              = 1/6 * dE / dt;
                trackInfo.lambda         = lambda;

                params.minStates       = minStates;
                params.maxStates       = maxStates;
                params.numFeatures     = numFeatures;
                params.numReinitialize = numReinitialize;
                params.numPerturbation = numPerturb;
                params.converged       = convergence;
                params.maxiter         = maxiter;
                params.verbose         = 1;

                deltaX = cell(trackInfo.numberOfTracks, 1);
                for i = 1:trackInfo.numberOfTracks
                    deltaX{i} = diff(splitX{i});
                end

                [trackInfo.vacf_exp, trackInfo.xbar_exp] = ...
                    CovarianceProperties(deltaX, numFeatures);

                tic;
                results = pEMv2_SPT(deltaX, trackInfo, params);
                toc;

                % Store inputs needed for bootstrap CI (computed separately via
                % computeBootstrapCI — not run here so a parallel pool failure
                % cannot prevent pEM results from being saved).
                if ~isSubset
                    obj.pEMBootstrapInputs.deltaX     = deltaX;
                    obj.pEMBootstrapInputs.rawResults = results;
                    obj.pEMBootstrapInputs.trackInfo  = trackInfo;
                    obj.pEMBootstrapInputs.params     = params;
                    obj.IsBootstrapComputed           = false;
                end

                vacf              = results.optimalVacf;
                results.optimalD  = (vacf(:,1) + 2*vacf(:,2)) / (2*trackInfo.dt);
                results.optimalS  = sqrt((vacf(:,1) - 2*results.optimalD*dt*(1-2*trackInfo.R)) / 2);

                [D, idx] = sort(results.optimalD);

                pEMTable = table;
                pEMTable.numRawTracks{1}   = numel(X);
                pEMTable.numSplitTracks(1) = results.trackInfo.numberOfTracks;
                pEMTable.trackInfo{1}      = results.trackInfo;
                pEMTable.params{1}         = results.params;
                pEMTable.optimalSize(1)    = results.optimalSize(1);
                pEMTable.optimalD{1}       = D;
                pEMTable.optimalS{1}       = results.optimalS(idx);
                pEMTable.optimalVacf{1}    = results.optimalVacf(idx, :);
                pEMTable.optimalP{1}       = results.optimalP(idx);
                pEMTable.optimalL(1)       = results.optimalL;
                pEMTable.BIC{1}            = results.BIC(idx);
                pEMTable.splitX{1}         = splitX;
                pEMTable.trackID{1}        = trackID;
                pEMTable.splitID{1}        = results.trackInfo.splitIndex;
                pEMTable.posteriorProb{1}  = results.posteriorProb(:, idx);

                posteriorProb              = pEMTable.posteriorProb{1};
                [maxPosteriorProb, optimalState] = max(posteriorProb, [], 2);

                vacf    = pEMTable.trackInfo{1}.vacf_exp;
                dt_tbl  = pEMTable.trackInfo{1}.dt;
                track_D = mean(vacf(:,1,:) + 2*vacf(:,2,:), 3) / 2 / dt_tbl;

                pEMTable.maxPosteriorProb{1} = maxPosteriorProb;
                pEMTable.optimalState{1}     = optimalState;
                pEMTable.track_D{1}          = track_D;

                maxPosteriorProb = pEMTable.maxPosteriorProb{1};
                optimalState     = pEMTable.optimalState{1};

                proportions = table;
                for i = 1:pEMTable.optimalSize
                    proportions.State(i)      = i;
                    proportions.Proportion(i) = length(find(optimalState==i)) / length(optimalState);
                    proportions.OptimalD(i)   = pEMTable.optimalD{1}(i);
                    proportions.OptimalS(i)   = pEMTable.optimalS{1}(i);
                end

                fprintf('Population fraction (global):\n');
                proportions

                if ~isSubset
                    obj.pEMParams.trackInfo    = trackInfo;
                    obj.pEMParams.params       = params;
                    obj.IspEMComputed          = true;
                    obj.pEMResults.pEMTable    = pEMTable;
                    obj.pEMResults.proportions = proportions;
                end
            end
            results = obj.pEMResults;
        end

        %% Bootstrap CI (separate from pEM)
        function ciResults = computeBootstrapCI(obj, varargin)
            % COMPUTEBOOTSTRAPCI  Compute confidence intervals on a completed pEM run.
            %
            %   ciResults = tc.computeBootstrapCI()
            %   ciResults = tc.computeBootstrapCI('nBoot', 200, 'parallel', false)
            %
            %   Must be called after getBayesianDiffusivity() has completed for
            %   the full collection (not a condition subset).  Results are stored
            %   in pEMResults.ciresults and included in subsequent calls to
            %   getBayesianDiffusivity().
            %
            %   Name-value options:
            %     'nBoot'          - bootstrap resamples          (default 100)
            %     'covModel'       - 'normal','confined','fBM'    (default 'normal')
            %     'alpha'          - CI level: 1-alpha            (default 0.05)
            %     'nRandomStarts'  - restarts per replicate        (default 5)
            %     'resampleLevel'  - 'parent' or 'subtracks'      (default 'parent')
            %     'parallel'       - use parfor                   (default true)
            %     'verbose'        - print progress               (default true)

            if ~obj.IspEMComputed
                error('TrajectoryCollection:pEMNotComputed', ...
                    'Run getBayesianDiffusivity() on the full collection first.');
            end
            if ~isfield(obj.pEMBootstrapInputs, 'deltaX') || isempty(obj.pEMBootstrapInputs.deltaX)
                error('TrajectoryCollection:noBootstrapInputs', ...
                    ['Bootstrap inputs not found. This can happen if pEM was run on ' ...
                     'a condition subset only. Run getBayesianDiffusivity() on the ' ...
                     'full collection (no Condition argument) first.']);
            end

            p = inputParser;
            addParameter(p, 'nBoot',         100,      @isnumeric);
            addParameter(p, 'covModel',      'normal',  @ischar);
            addParameter(p, 'alpha',          0.05,     @isnumeric);
            addParameter(p, 'nRandomStarts',  5,        @isnumeric);
            addParameter(p, 'resampleLevel', 'parent',  @ischar);
            addParameter(p, 'parallel',       true,     @islogical);
            addParameter(p, 'verbose',        true,     @islogical);
            parse(p, varargin{:});
            o = p.Results;

            deltaX     = obj.pEMBootstrapInputs.deltaX;
            rawResults = obj.pEMBootstrapInputs.rawResults;
            trackInfo  = obj.pEMBootstrapInputs.trackInfo;
            params     = obj.pEMBootstrapInputs.params;

            fprintf('\n=== Bootstrap CI ===\n');
            fprintf('Timestamp  : %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
            fprintf('nBoot      : %d\n', o.nBoot);
            fprintf('covModel   : %s\n', o.covModel);
            fprintf('alpha      : %.3f (%.0f%% CI)\n', o.alpha, 100*(1-o.alpha));
            fprintf('resample   : %s\n', o.resampleLevel);
            fprintf('parallel   : %s\n', ifelse(o.parallel, 'yes', 'no'));

            tic;
            ciResults = pEMv2_bootstrap_CI(deltaX, rawResults, trackInfo, params, ...
                'nBoot',         o.nBoot, ...
                'covModel',      o.covModel, ...
                'alpha',         o.alpha, ...
                'nRandomStarts', o.nRandomStarts, ...
                'resampleLevel', o.resampleLevel, ...
                'parallel',      o.parallel, ...
                'verbose',       o.verbose);
            toc;

            obj.pEMResults.ciresults = ciResults;
            obj.IsBootstrapComputed  = true;
            fprintf('Bootstrap CI stored in pEMResults.ciresults.\n');
            fprintf('===================\n\n');
        end

        %% Lifetime
        function results = getLifetime(obj, varargin)
            if ~obj.IsLifetimeComputed
                p = inputParser;
                addParameter(p, 'rmin',           0.25);
                addParameter(p, 'rmax',           0.4);
                addParameter(p, 'minTrackLength', 4);
                addParameter(p, 'maxTrackLength', 400);
                addParameter(p, 'UseRaw',         false);
                parse(p, varargin{:});

                if p.Results.UseRaw
                    tracks = obj.getAllRawTracks();
                else
                    tracks = obj.getAllCulledTracks();
                end
                if isempty(tracks)
                    error('No tracks available for lifetime computation.');
                end

                tracklen = cellfun('size', tracks, 1);
                X = cell(1, numel(tracks));
                for itrack = 1:numel(tracks)
                    X{itrack} = tracks{itrack}(:, 1:2);   % already in µm
                end

                idx  = find(tracklen > p.Results.minTrackLength & ...
                            tracklen < p.Results.maxTrackLength);
                tmp{1} = X(idx);
                [tracklen, x, f, flo, fup, jump_hist] = ...
                    find_bound_particles_2(tmp, p.Results.rmin, p.Results.rmax);
                prctile(jump_hist, 99)

                obj.LifetimeResults.x           = x;
                obj.LifetimeResults.f           = f;
                obj.LifetimeResults.flo         = flo;
                obj.LifetimeResults.fup         = fup;
                obj.LifetimeResults.TrackLengths = tracklen;
                obj.IsLifetimeComputed = true;
            end
            results = obj.LifetimeResults;
        end

        function obj = recomputeLifetime(obj, varargin)
            obj.IsLifetimeComputed = false;
            obj.LifetimeResults    = struct();
            obj.getLifetime(varargin{:});
            fprintf('Lifetime recomputed.\n');
        end

        function obj = recomputeRL(obj, varargin)
            obj.IsRLComputed = false;
            obj.RLResults    = struct();
            obj.getRLDecomposition(varargin{:});
            fprintf('RL decomposition recomputed.\n');
        end

        function obj = recomputeMSD(obj, varargin)
            obj.IsMSDComputed = false;
            obj.MSDResults    = struct();
            obj.getMSD(varargin{:});
            fprintf('MSD recomputed.\n');
        end

        %% Visualization
        function fig = visualizeSubsets(obj, varargin)
            p = inputParser;
            addParameter(p, 'Condition', '');
            addParameter(p, 'FileID',    '');
            addParameter(p, 'NumTracks', inf);
            addParameter(p, 'UseRaw',    false);
            parse(p, varargin{:});

            idx = true(size(obj.Wrappers));
            if ~isempty(p.Results.Condition)
                idx = idx & strcmp(obj.Metadata.Condition, p.Results.Condition);
            end
            if ~isempty(p.Results.FileID)
                idx = idx & strcmp(obj.Metadata.FileID, p.Results.FileID);
            end
            subWrappers = obj.Wrappers(idx);

            fig    = figure;
            hold on;
            colors = lines(length(subWrappers));
            for w = 1:length(subWrappers)
                tracks = ifelse(p.Results.UseRaw, subWrappers{w}.getRawTracks(), ...
                    subWrappers{w}.getCulledTracks());
                if length(tracks) > p.Results.NumTracks
                    tracks = tracks(randperm(length(tracks), p.Results.NumTracks));
                end
                for t = 1:length(tracks)
                    plot(tracks{t}(:,1), tracks{t}(:,2), 'Color', colors(w,:));
                end
                if ~isempty(subWrappers{w}.getNucleusMask())
                    contour(subWrappers{w}.getNucleusMask(), [0.5 0.5], 'k--');
                end
            end
            title('Trajectory Subsets');
            xlabel('X'); ylabel('Y');
            hold off;
        end

        function figs = plotMSD(obj, varargin)
            % PLOTMSD  Plot ensemble MSD with bootstrap CI and median fBM fit.
            %
            %   figs = tc.plotMSD()
            %   figs = tc.plotMSD('Color',  [0 0.4 0.8])
            %   figs = tc.plotMSD('Title',  'H2B Control')
            %   figs = tc.plotMSD('ShowBootstraps', true)   % overlay individual fits
            %
            % Returns a figure handle.  Call tc.getMSD() first.

            if ~obj.IsMSDComputed
                error('plotMSD: run tc.getMSD() first.');
            end

            p = inputParser;
            addParameter(p, 'Color',          [0 0.4 0.8], @isnumeric);
            addParameter(p, 'Title',          '',          @(x) ischar(x)||isstring(x));
            addParameter(p, 'ShowBootstraps', false,       @islogical);
            parse(p, varargin{:});

            r      = obj.MSDResults;
            tau    = r.lag_axis;
            emsd   = r.ensemble_mean;
            bmeans = r.boot_means;
            bfits  = r.boot_fits;
            dt     = r.dt;
            frac   = r.frac;
            lags   = r.lag_frames;
            clr    = p.Results.Color;

            ci_lo = prctile(bmeans, 2.5,  2);
            ci_hi = prctile(bmeans, 97.5, 2);

            figs = figure('Name', 'Ensemble MSD');
            patch([tau; flipud(tau)], [ci_lo; flipud(ci_hi)], clr, ...
                  'EdgeColor', 'none', 'FaceAlpha', 0.25);
            hold on;

            % Individual bootstrap fit curves (optional)
            if p.Results.ShowBootstraps
                for b = 1:size(bfits, 1)
                    fp = bfits(b, :);
                    if any(isnan(fp)), continue; end
                    fc = msd_model_curve(fp, lags, dt, frac);
                    plot(tau, fc, 'Color', [clr 0.08], 'LineWidth', 0.5);
                end
            end

            % Ensemble mean
            plot(tau, emsd, '-', 'Color', clr, 'LineWidth', 2.5);

            % Median fit curve
            valid = bfits(~any(isnan(bfits), 2), :);
            if ~isempty(valid)
                med_fp   = median(valid, 1);
                fit_crv  = msd_model_curve(med_fp, lags, dt, frac);
                plot(tau, fit_crv, '--', 'Color', clr * 0.6, 'LineWidth', 1.5);
                med_alpha = med_fp(3);
                ci_alpha  = prctile(valid(:,3), [2.5 97.5]);
                legend({'95% CI', 'Ensemble mean', ...
                        sprintf('fBM fit  \\alpha=%.2f [%.2f–%.2f]', ...
                            med_alpha, ci_alpha(1), ci_alpha(2))}, ...
                       'Location', 'northwest', 'Box', 'off');
            else
                legend({'95% CI', 'Ensemble mean'}, 'Location', 'northwest', 'Box', 'off');
            end

            set(gca, 'XScale', 'log', 'YScale', 'log', ...
                     'FontSize', 18, 'FontWeight', 'bold', 'LineWidth', 1.5);
            xlabel('\tau (s)');
            ylabel('MSD (\mum^2)');
            ttl = char(p.Results.Title);
            if isempty(ttl)
                ttl = sprintf('Ensemble MSD  (n=%d tracks)', r.n_tracks);
            end
            title(ttl);
            axis square; box off;
        end

        function figs = plotRLDecomposition(obj, varargin)
            % PLOTRLDECOMPOSITION  Reproduce the standard RL analysis figures
            %                      from stored results without re-running RL_HILO.
            %
            %   figs = tc.plotRLDecomposition()
            %   figs = tc.plotRLDecomposition('Title', 'H2B Control')
            %
            % Produces three figures:
            %   Fig 1 – van Hove correlation G(r,τ): data + RL fit
            %   Fig 2 – P(MSD): Richardson-Lucy deconvolved MSD distribution
            %   Fig 3 – Per-state ensemble MSD with 99% CI error bars
            %
            % Call tc.getRLDecomposition() first.

            if ~obj.IsRLComputed
                error('plotRLDecomposition: run tc.getRLDecomposition() first.');
            end

            p = inputParser;
            addParameter(p, 'Title', '', @(x) ischar(x)||isstring(x));
            parse(p, varargin{:});

            cq      = obj.RLResults.computed_quantities;
            ttl     = char(p.Results.Title);
            if isempty(ttl), ttl = cq.ident; end

            dt      = obj.getFrameIntervalForCondition([]);
            n_lags  = size(cq.msd, 1);
            tau_rl  = (1:n_lags)' * dt;
            n_states = numel(cq.classified_tracks);
            colors   = lines(n_states);

            % --- Fig 1: van Hove correlation ---------------------------------
            figs(1) = figure('Name', 'van Hove correlation');
            plot(cq.x, cq.vanHove, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'none');
            hold on;
            plot(cq.x, cq.Gs, 'r-', 'LineWidth', 2);
            set(gca, 'YScale', 'log', 'FontSize', 24, 'FontWeight', 'bold', 'LineWidth', 1);
            xlabel('r  (\mum)');
            ylabel('G(r,\tau)');
            axis([0 1 1e-4 100]);
            axis square; box off;
            title(ttl);
            legend({'Data', 'RL fit'}, 'Location', 'northeast', 'Box', 'off');

            % --- Fig 2: P(MSD) -----------------------------------------------
            figs(2) = figure('Name', 'P(MSD)');
            semilogx(obj.RLResults.M, cq.P1norm, 'k-', 'LineWidth', 2);
            set(gca, 'FontSize', 24, 'FontWeight', 'bold', 'LineWidth', 1);
            xlabel('MSD  (\mum^2)');
            ylabel('P(MSD)');
            axis([1e-3 2e-1 0 max(cq.P1norm) * 1.1]);
            title(ttl);
            box off;

            % --- Fig 3: per-state ensemble MSD -------------------------------
            figs(3) = figure('Name', 'Per-state MSD');
            hold on;
            for s = 1:n_states
                n_s = numel(cq.classified_tracks{s});
                if n_s == 0, continue; end
                errorbar(tau_rl, cq.msd(:, s), cq.msderr(:, s), ...
                         '-o', 'Color', colors(s, :), 'LineWidth', 1.5, ...
                         'MarkerSize', 4, 'CapSize', 3, ...
                         'DisplayName', sprintf('State %d  (n=%d)', s, n_s));
            end
            set(gca, 'XScale', 'log', 'YScale', 'log', ...
                     'FontSize', 18, 'FontWeight', 'bold', 'LineWidth', 1.5);
            xlabel('\tau (s)');
            ylabel('MSD  (\mum^2)');
            title([ttl '  —  per-state MSD']);
            legend('show', 'Location', 'northwest', 'Box', 'off');
            axis square; box off;
        end

        function globalfig = plotpEMstats(obj)
            if ~obj.IspEMComputed
                warning('First compute pEM');
                return;
            end
            proportions  = obj.pEMResults.proportions;
            pp           = obj.pEMResults.pEMTable.posteriorProb{1};
            max_pp       = max(pp, [], 2);
            sort_pp      = sort(pp, 2, 'descend');
            delta_pp     = sort_pp(:,1) - sort_pp(:,2);
            idx          = find(max_pp > 0.5 & delta_pp > 0.2);
            optimalState = obj.pEMResults.pEMTable.optimalState{1};
            optimalState(setdiff(1:length(optimalState), idx)) = 9;
            actual_num   = length(optimalState) - length(find(optimalState == 9));
            kept_states  = true(height(proportions), 1);
            for ii = 1:height(proportions)
                if length(find(optimalState==ii)) / actual_num < 0.025
                    optimalState(find(optimalState==ii)) = 9;
                    kept_states(ii) = false;
                else
                    sprintf('state %d proportion %f', ii, length(find(optimalState==ii)) / actual_num);
                end
            end
            splitLength       = 7;
            colors            = [[0.7 0.3 0.3];[1 0 0];[0 0 1];[0 0.8 0.2];[0 1 0];[0.2 0.2 0.2]];
            msd_all_states    = zeros(splitLength-1, height(proportions));
            var_all_states    = zeros(splitLength-1, height(proportions));
            for ii = 1:length(kept_states)
                if kept_states(ii)
                    idx2 = find(optimalState == ii);
                    dum  = obj.pEMResults.pEMTable.splitX{1}(idx2);
                    ensemble_msd  = zeros(splitLength-1, length(idx2));
                    ensemble_lags = zeros(splitLength-1, length(idx2));
                    for jj = 1:length(dum)
                        [msd, num_lags_tracks, ~] = simple_msd(dum{jj}(:,1:2), (1:7)');
                        ensemble_msd(:,jj)  = msd;
                        ensemble_lags(:,jj) = num_lags_tracks;
                    end
                    msd_all_states(:,ii) = mean(ensemble_msd ./ ensemble_lags, 2);
                    var_all_states(:,ii) = std(ensemble_msd ./ ensemble_lags, 0, 2) ./ ...
                                           sqrt(sum(ensemble_lags, 2));
                end
            end
            globalfig(1) = figure('Position', [1 1 0.5 1] .* get(0,'Screensize'), 'Visible','on');
            dt = obj.getFrameIntervalForCondition([]);
            for ii = 1:length(kept_states)
                if kept_states(ii)
                    errorbar(dt*(0:6), [0; msd_all_states(:,ii)], [0; var_all_states(:,ii)], ...
                             's', 'Color', colors(ii,:), 'markersize', 8);
                    hold on;
                    plot(dt*(0:6), [0; msd_all_states(:,ii)], 'Color', colors(ii,:), 'linewidth', 1);
                end
            end
            axis square; box off;
            set(gca, 'linewidth', 2, 'fontsize', 24, 'fontweight', 'bold');
        end

        function globalfig = plotLifetime(obj)
            x   = obj.LifetimeResults.x;
            f   = obj.LifetimeResults.f;
            flo = obj.LifetimeResults.flo;
            fup = obj.LifetimeResults.fup;
        end

        %% Metadata / query helpers
        function obj = setCondition(obj, fileID, newCondition)
            idx = strcmp(obj.Metadata.FileID, fileID);
            if any(idx)
                obj.Metadata.Condition(idx) = newCondition;
            else
                warning('FileID not found.');
            end
        end

        function tracklen = getTrackLengthDistribution(obj)
            tracks   = obj.getAllCulledTracks();
            tracklen = cellfun('size', tracks, 1);
        end

        function tracks = getTracksByCondition(obj, condition, varargin)
            % GETTRACKSBYCONDITION  Return culled (or raw) tracks for a condition.
            %
            %   tracks = tc.getTracksByCondition('Control')
            %   tracks = tc.getTracksByCondition('Treated', 'UseRaw', true)

            p = inputParser;
            addRequired(p, 'condition');
            addParameter(p, 'UseRaw', false, @islogical);
            parse(p, condition, varargin{:});

            if ischar(condition) || isstring(condition)
                condition = {char(condition)};
            end

            matchIdx = false(size(obj.Metadata.Condition));
            for c = condition(:)'
                matchIdx = matchIdx | strcmp(obj.Metadata.Condition, c);
            end

            if ~any(matchIdx)
                warning('No files found for condition(s): %s', strjoin(condition, ', '));
                tracks = {};
                return;
            end

            matchingWrappers = obj.Wrappers(matchIdx);
            tracks = cell(0,1);

            for i = 1:length(matchingWrappers)
                tw = matchingWrappers{i};
                if p.Results.UseRaw
                    thisTracks = tw.getRawTracks();
                else
                    thisTracks = tw.getCulledTracks();
                end
                if ~isempty(thisTracks)
                    if ~iscolumn(thisTracks), thisTracks = thisTracks(:); end
                    tracks = [tracks; thisTracks]; %#ok<AGROW>
                end
            end

            fprintf('Returned %d tracks for condition(s): %s\n', ...
                length(tracks), strjoin(condition, ', '));
        end

        function tracks = getRelativeTracksByCondition(obj, condition)
            % GETRELATIVETRACKSBYCONDITION  Return relative tracks for a condition.
            if ischar(condition) || isstring(condition)
                condition = {char(condition)};
            end
            matchIdx = false(size(obj.Metadata.Condition));
            for c = condition(:)'
                matchIdx = matchIdx | strcmp(obj.Metadata.Condition, c);
            end
            if ~any(matchIdx)
                warning('No files found for condition(s): %s', strjoin(condition, ', '));
                tracks = {};
                return;
            end
            tracks = cell(0, 1);
            for tw = obj.Wrappers(matchIdx)
                rt = tw{1}.getRelativeTracks();
                if ~isempty(rt)
                    if ~iscolumn(rt), rt = rt(:); end
                    tracks = [tracks; rt]; %#ok<AGROW>
                end
            end
        end

        function summary(obj)
            fprintf('\n=== TrajectoryCollection Summary ===\n');
            fprintf('Number of files/wrappers: %d\n', length(obj.Wrappers));

            if ~isempty(obj.Wrappers)
                rawCounts    = cellfun(@(tw) tw.getNumRawTracks(),    obj.Wrappers);
                culledCounts = cellfun(@(tw) tw.getNumCulledTracks(), obj.Wrappers);
                relCounts    = cellfun(@(tw) numel(tw.getRelativeTracks()), obj.Wrappers);
                srcTypes     = cellfun(@(tw) tw.getSourceType(),       obj.Wrappers, 'UniformOutput', false);
                fis          = cellfun(@(tw) tw.getFrameInterval(),    obj.Wrappers);

                fprintf('Total raw tracks     : %d\n', sum(rawCounts));
                fprintf('Total culled tracks  : %d\n', sum(culledCounts));
                fprintf('Total relative tracks: %d\n', sum(relCounts));

                fprintf('\nPer-wrapper details:\n');
                for w = 1:length(obj.Wrappers)
                    fiStr = ifelse(isnan(fis(w)), 'NOT SET', sprintf('%.4f s', fis(w)));
                    fprintf('  [%s] source=%-4s  culled=%4d  px=%.4f  dt=%s\n', ...
                        char(obj.Metadata.FileID(w)), srcTypes{w}, ...
                        culledCounts(w), obj.Wrappers{w}.getPixelSize(), fiStr);
                end
            else
                fprintf('Total raw tracks     : 0\n');
                fprintf('Total culled tracks  : 0\n');
                fprintf('Total relative tracks: 0\n');
            end

            uniqueConditions = unique(obj.Metadata.Condition);
            if ~isempty(uniqueConditions) && ~all(strcmp(uniqueConditions, ''))
                fprintf('Conditions           : %s\n', strjoin(uniqueConditions, ', '));
            else
                fprintf('Conditions           : (none specified)\n');
            end

            fprintf('RL decomp computed   : %s\n', ifelse(obj.IsRLComputed,        'Yes', 'No'));
            fprintf('MSD computed         : %s\n', ifelse(obj.IsMSDComputed,       'Yes', 'No'));
            fprintf('Bayesian computed    : %s\n', ifelse(obj.IspEMComputed,       'Yes', 'No'));
            fprintf('Bootstrap CI computed: %s\n', ifelse(obj.IsBootstrapComputed, 'Yes', 'No'));
            fprintf('====================================\n\n');
        end

        function exportToCSV(obj, outputFile, varargin)
            p = inputParser;
            addParameter(p, 'UseRaw', false);
            parse(p, varargin{:});

            tracks  = ifelse(p.Results.UseRaw, obj.getAllRawTracks(), obj.getAllCulledTracks());
            allData = [];
            for i = 1:length(tracks)
                meta    = repmat(obj.Metadata(i,:), size(tracks{i},1), 1);
                allData = [allData; [tracks{i}, meta.FileID, meta.Condition]]; %#ok<AGROW>
            end
            writetable(array2table(allData), outputFile);
        end

        function tw = getWrapperForTrack(obj, trackIdx, varargin)
            p = inputParser;
            addParameter(p, 'UseRaw', false, @islogical);
            parse(p, varargin{:});

            tracks     = ifelse(p.Results.UseRaw, obj.getAllRawTracks(), obj.getAllCulledTracks());
            cumulative = cumsum(cellfun(@(w) w.getNumCulledTracks(), obj.Wrappers));
            if p.Results.UseRaw
                cumulative = cumsum(cellfun(@(w) w.getNumRawTracks(), obj.Wrappers));
            end

            wrapperIdx = find(cumulative >= trackIdx, 1, 'first');
            if isempty(wrapperIdx)
                error('Track index %d out of range (max %d).', trackIdx, length(tracks));
            end
            tw = obj.Wrappers{wrapperIdx};
        end

        function singleTrack = getTrack(obj, trackIdx, varargin)
            p = inputParser;
            addParameter(p, 'UseRaw', false);
            parse(p, varargin{:});

            allTracks = ifelse(p.Results.UseRaw, obj.getAllRawTracks(), obj.getAllCulledTracks());
            if trackIdx < 1 || trackIdx > length(allTracks)
                error('Track index %d out of range (1 to %d).', trackIdx, length(allTracks));
            end
            singleTrack = allTracks{trackIdx};
        end

        %% Save
        function save(obj, filename)
            % SAVE  Persist the TrajectoryCollection to a .mat file.
            %
            %   tc.save('results.mat')
            %
            %   Reload with:
            %     data = load('results.mat');  tc = data.tc;
            %     tc.reloadTracks()            % optional: rebuild track data from CSVs
            %
            %   Track matrices (RawTracks / RelativeTracks) are automatically
            %   stripped from each Wrapper at serialization time via saveobj —
            %   the live in-memory object is NOT modified.  All analysis results
            %   (pEM, bootstrap CI, MSD, lifetime) and enough metadata to
            %   reconstruct tracks (file paths, ROIs, Parameters, PixelSize,
            %   FrameInterval) are preserved.
            %
            %   Typical file sizes:
            %     Without tracks : ~100-400 MB  (dominated by pEMBootstrapInputs)
            %     tc.reloadTracks() : 5-10 min to rebuild from original CSVs

            if ~endsWith(filename, '.mat')
                filename = [filename '.mat'];
            end

            % Clear TC-level track caches — they are derived from Wrapper tracks
            % and will be stale after a save/load cycle anyway.
            obj.AllRawTracks           = {};
            obj.AllCulledTracks        = {};
            obj.AllRelativeTracks      = {};
            obj.IsAllRawCollected      = false;
            obj.IsAllCulledCollected   = false;
            obj.IsAllRelativeCollected = false;

            tc = obj; %#ok<NASGU>
            saveTCToFile(filename, tc);
            fprintf('Saved TrajectoryCollection (%d wrappers) to %s\n', ...
                numel(obj.Wrappers), filename);
        end

        %% Reload tracks from original CSV files
        function reloadTracks(obj, varargin)
            % RELOADTRACKS  Rebuild RawTracks in each Wrapper from original CSV files.
            %
            %   tc.reloadTracks()
            %   tc.reloadTracks('OldRoot', '/old/mount/', 'NewRoot', '/new/mount/')
            %
            %   Call this after loading a saved TrajectoryCollection to restore
            %   track data.  Analysis results (pEM, MSD, etc.) are unaffected.
            %
            %   If CSV files were added with absolute paths (recommended), no
            %   arguments are needed — reloadTracks() finds them automatically.
            %
            %   Name-value options:
            %     'OldRoot' / 'NewRoot' - pair that replaces a path prefix when
            %                    data has moved (e.g. different mount point).
            %                    Both must be supplied together.
            %                    E.g.: 'OldRoot','/lab/server/', 'NewRoot','/Volumes/Lab/'
            %                    preserves subdirectory structure under the new root.

            p = inputParser;
            addParameter(p, 'OldRoot', '', @ischar);
            addParameter(p, 'NewRoot', '', @ischar);
            parse(p, varargin{:});
            oldRoot = p.Results.OldRoot;
            newRoot = p.Results.NewRoot;

            nReloaded = 0;
            for i = 1:numel(obj.Wrappers)
                tw = obj.Wrappers{i};

                if ~isempty(tw.getRawTracks())
                    continue;   % already has tracks
                end

                if ~strcmp(tw.SourceType, 'csv')
                    warning('TrajectoryCollection:reloadTracks', ...
                        'Wrapper %d is SMD-sourced; tracks cannot be reloaded from CSV.', i);
                    continue;
                end

                fname = tw.FileName;
                if ~isempty(oldRoot) && ~isempty(newRoot) && startsWith(fname, oldRoot)
                    fname = fullfile(newRoot, fname(numel(oldRoot)+1:end));
                end

                if isempty(fname) || ~isfile(fname)
                    warning('TrajectoryCollection:reloadTracks', ...
                        'File not found for wrapper %d: %s', i, fname);
                    continue;
                end

                wasRelative = tw.IsRelative;
                tw.readData(fname);
                % ROIs were preserved by saveobj/loadobj — no need to re-run loadNucleusMask
                if tw.IsCulled || ~isempty(tw.ROIs)
                    tw.cull_tracks();
                    if wasRelative
                        tw.remove_relative_motion();
                    end
                end
                nReloaded = nReloaded + 1;
            end

            obj.invalidateCollections();
            fprintf('Reloaded tracks for %d/%d wrappers.\n', nReloaded, numel(obj.Wrappers));
        end

    end % public methods

    methods (Access = private)

        function collectAllTracks(obj)
            nWrappers = length(obj.Wrappers);
            if nWrappers == 0
                obj.AllRawTracks       = {};
                obj.AllCulledTracks    = {};
                obj.AllRelativeTracks  = {};
                obj.IsAllRawCollected      = true;
                obj.IsAllCulledCollected   = true;
                obj.IsAllRelativeCollected = true;
                return;
            end

            raw      = cell(0,1);
            culled   = cell(0,1);
            relative = cell(0,1);

            for i = 1:nWrappers
                tw           = obj.Wrappers{i};
                thisRaw      = tw.getRawTracks();
                thisCulled   = tw.getCulledTracks();
                thisRelative = tw.getRelativeTracks();

                if ~isempty(thisRaw)      && ~iscolumn(thisRaw),      thisRaw      = thisRaw(:);      end
                if ~isempty(thisCulled)   && ~iscolumn(thisCulled),   thisCulled   = thisCulled(:);   end
                if ~isempty(thisRelative) && ~iscolumn(thisRelative),  thisRelative = thisRelative(:); end

                raw      = [raw;      thisRaw];      %#ok<AGROW>
                culled   = [culled;   thisCulled];   %#ok<AGROW>
                relative = [relative; thisRelative]; %#ok<AGROW>
            end

            obj.AllRawTracks       = raw;
            obj.AllCulledTracks    = culled;
            obj.AllRelativeTracks  = relative;
            obj.IsAllRawCollected      = true;
            obj.IsAllCulledCollected   = true;
            obj.IsAllRelativeCollected = true;

            fprintf('Collected %d raw, %d culled, and %d relative tracks across %d FOVs.\n', ...
                length(raw), length(culled), length(relative), nWrappers);
        end

        function dt = getFrameIntervalForCondition(obj, condition)
            % Returns the frame interval (seconds) for the requested condition.
            % Errors if any matching wrapper has FrameInterval = NaN.
            % Warns if wrappers disagree on the value.

            if isempty(condition)
                wrappers = obj.Wrappers;
            else
                if ischar(condition) || isstring(condition)
                    condition = {char(condition)};
                end
                matchIdx = false(size(obj.Metadata.Condition));
                for c = condition(:)'
                    matchIdx = matchIdx | strcmp(obj.Metadata.Condition, c);
                end
                wrappers = obj.Wrappers(matchIdx);
            end

            if isempty(wrappers)
                error('TrajectoryCollection:noWrappers', ...
                      'No wrappers found for the specified condition.');
            end

            fis = cellfun(@(tw) tw.getFrameInterval(), wrappers);

            if any(isnan(fis))
                error('TrajectoryCollection:frameIntervalNotSet', ...
                    ['FrameInterval is NaN for %d wrapper(s). ' ...
                     'Call tw.readParams(tomlFile) or set FrameInterval explicitly ' ...
                     'before running analysis.'], sum(isnan(fis)));
            end

            uniqueFIs = unique(fis);
            if numel(uniqueFIs) > 1
                warning('TrajectoryCollection:inconsistentFrameInterval', ...
                    'Wrappers have %d different FrameInterval values: [%s]. Using median = %.4f s.', ...
                    numel(uniqueFIs), num2str(uniqueFIs(:)', '%.4f '), median(fis));
            end

            dt = median(fis);
        end

        function invalidateCollections(obj)
            obj.IsAllRawCollected      = false;
            obj.IsAllCulledCollected   = false;
            obj.IsAllRelativeCollected = false;
            obj.IsRLComputed           = false;
            obj.IsMSDComputed          = false;
            obj.IspEMComputed          = false;
            obj.IsBootstrapComputed    = false;
            obj.AllRawTracks       = {};
            obj.AllCulledTracks    = {};
            obj.AllRelativeTracks  = {};
            obj.RLResults          = struct();
            obj.MSDResults         = struct();
            obj.pEMResults             = struct();
            obj.pEMBootstrapInputs     = struct();
        end

        function addGlobalParameters(obj, params)
            for w = obj.Wrappers
                w{1}.Parameters = mergeStructs(w{1}.Parameters, params);
            end
        end

    end % private methods
end

% -------------------------------------------------------------------------
% Module-level helpers
% -------------------------------------------------------------------------
function str = ifelse(cond, t, f)
    if cond, str = t; else, str = f; end
end

function curve = msd_model_curve(fp, lags, dt, frac)
% Reconstruct fBM MSD model values from fit parameters.
%   fp = [G, sigma_sq, alpha]; lags = integer lag indices.
G = fp(1); sig2 = fp(2); alpha = fp(3);
b     = (abs(1 + frac./lags).^(2+alpha) + abs(1-frac./lags).^(2+alpha) - 2) ./ (frac./lags).^2;
curve = G / ((1+alpha)*(2+alpha)) .* ((dt.*lags).^alpha .* b - 2*(frac*dt)^alpha) + 2*sig2;
end

function s = mergeStructs(s1, s2)
    s      = s1;
    fields = fieldnames(s2);
    for i = 1:length(fields)
        s.(fields{i}) = s2.(fields{i});
    end
end
