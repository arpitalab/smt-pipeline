function results = classifyTracksGMM(tracks, varargin)
%CLASSIFYTRACKSGMM  Classify tracks via 2-D GMM in (log10(Pc), alpha) space.
%
%   results = classifyTracksGMM(tracks)
%   results = classifyTracksGMM(tracks, 'dt', 0.2, 'SigmaLoc', 0.040, ...)
%
%   Computes per-track packing fraction (Pc) and per-track fBM anomalous
%   exponent (α) for each qualifying track, then fits a Gaussian mixture
%   model (GMM) in the 2-D feature space [log10(Pc), α].  The number of
%   components is chosen by BIC over the range MinComponents:MaxComponents.
%
%   Inputs:
%     tracks  - Cell array of Nx2 (or Nx3) trajectory matrices (µm).
%               Rows = frames, col 1 = x, col 2 = y.
%
%   Name-value options (packing fraction):
%     'WindowSize'     - Sliding window for Pc (frames, default: 20).
%     'MinTrackLength' - Skip tracks shorter than this (default: WindowSize+1).
%     'MinArea'        - Min convex-hull area to accept a window (µm², default: 1e-4).
%
%   Name-value options (fBM per-track fit):
%     'dt'             - Frame interval (s, default: 0.2).
%     'ExposureFraction' - Te/dt for fBM kernel (default: 1).
%     'SubtrackLength' - Subtrack length for fBM population MLE (default: 20).
%     'SigmaLoc'       - Localization precision (µm, default: 0.040).
%                        Used as lower bound / prior; passed to fitFBM_pertracks.
%     'MinAlpha'       - Discard per-track fits with α < this (default: 0.05).
%     'MaxAlpha'       - Discard per-track fits with α > this (default: 2.0).
%
%   Name-value options (GMM):
%     'MinComponents'  - Minimum number of GMM components to try (default: 2).
%     'MaxComponents'  - Maximum number of GMM components to try (default: 4).
%     'Replicates'     - GMM fitting replicates for robustness (default: 10).
%     'CovType'        - GMM covariance type: 'full'|'diagonal' (default: 'full').
%
%   Name-value options (output):
%     'Plot'           - true/false — show scatter + GMM ellipses (default: true).
%     'Title'          - String for figure title (default: '').
%
%   Output struct fields:
%     .track_index     Qx1  indices into input tracks for qualifying tracks
%     .log10Pc         Qx1  log10 of per-track median packing fraction
%     .alpha           Qx1  per-track fBM α
%     .K               Qx1  per-track fBM K
%     .labels          Qx1  GMM component assignment (1..nComp)
%     .gmm             fitted gmdistribution object
%     .n_components    scalar  number of components chosen by BIC
%     .bic             vector  BIC values for MinComponents:MaxComponents
%     .classified_tracks  nComp x 1 cell of track index vectors (per state)

p = inputParser;
addRequired(p,  'tracks',                           @iscell);
% Packing fraction
addParameter(p, 'WindowSize',       20,             @isnumeric);
addParameter(p, 'MinTrackLength',   [],             @isnumeric);
addParameter(p, 'MinArea',          1e-4,           @isnumeric);
% fBM
addParameter(p, 'dt',               0.2,            @isnumeric);
addParameter(p, 'ExposureFraction', 1,              @isnumeric);
addParameter(p, 'SubtrackLength',   20,             @isnumeric);
addParameter(p, 'SigmaLoc',         0.040,          @isnumeric);
addParameter(p, 'MinAlpha',         0.05,           @isnumeric);
addParameter(p, 'MaxAlpha',         2.0,            @isnumeric);
% GMM
addParameter(p, 'MinComponents',    2,              @isnumeric);
addParameter(p, 'MaxComponents',    4,              @isnumeric);
addParameter(p, 'Replicates',       10,             @isnumeric);
addParameter(p, 'CovType',          'full',         @ischar);
% Output
addParameter(p, 'Plot',             true,           @islogical);
addParameter(p, 'Title',            '',             @ischar);
parse(p, tracks, varargin{:});
o = p.Results;

% ── 1. Per-track packing fraction ────────────────────────────────────────
fprintf('classifyTracksGMM: computing packing fraction...\n');
pc_res = computePackingFraction(tracks, ...
    'WindowSize',    o.WindowSize, ...
    'MinTrackLength', o.MinTrackLength, ...
    'MinArea',        o.MinArea, ...
    'Aggregate',      'median');

% pc_res.track_index contains indices into tracks of qualifying tracks
pc_idx    = pc_res.track_index;   % Px1
pc_vals   = pc_res.pc;            % Px1 median Pc per track

% ── 2. Per-track fBM α ───────────────────────────────────────────────────
fprintf('classifyTracksGMM: fitting per-track fBM α...\n');
% fitFBM_pertracks returns results indexed through qualifying tracks
fbm_res = fitFBM_pertracks(tracks, o.dt, ...
    'ExposureFraction',   o.ExposureFraction, ...
    'MaxSubtrackLength',  o.SubtrackLength, ...
    'MinTrackLength',     o.WindowSize + 1, ...
    'MinStepVar',         0, ...
    'MaxAlpha',           o.MaxAlpha, ...
    'Verbose',            false);

% Apply MinAlpha post-fit filter
alpha_keep = fbm_res.alpha >= o.MinAlpha;
fields_to_filter = fieldnames(fbm_res);
for fi = 1:numel(fields_to_filter)
    v = fbm_res.(fields_to_filter{fi});
    if isnumeric(v) && numel(v) == numel(alpha_keep)
        fbm_res.(fields_to_filter{fi}) = v(alpha_keep);
    end
end

% fbm_res.original_index: index into input tracks for each fitted track
fbm_idx   = fbm_res.original_index;  % Fx1
fbm_alpha = fbm_res.alpha;           % Fx1
fbm_K     = fbm_res.K;               % Fx1

% ── 3. Intersect: tracks that passed BOTH filters ─────────────────────────
[common_track_idx, ia, ib] = intersect(pc_idx, fbm_idx);
if isempty(common_track_idx)
    error('classifyTracksGMM: no tracks survived both Pc and fBM filters.');
end

log10Pc = log10(pc_vals(ia));
alpha   = fbm_alpha(ib);
K       = fbm_K(ib);

fprintf('classifyTracksGMM: %d tracks qualify for GMM (Pc∩fBM intersection).\n', ...
    numel(common_track_idx));

% Remove any NaN/Inf rows
valid = isfinite(log10Pc) & isfinite(alpha);
if any(~valid)
    warning('classifyTracksGMM: dropping %d tracks with NaN/Inf features.', sum(~valid));
    common_track_idx = common_track_idx(valid);
    log10Pc          = log10Pc(valid);
    alpha            = alpha(valid);
    K                = K(valid);
end

X = [log10Pc, alpha];   % N x 2 feature matrix

% ── 4. GMM with BIC model selection ──────────────────────────────────────
k_range = o.MinComponents : o.MaxComponents;
bic_vals = nan(1, numel(k_range));

fprintf('classifyTracksGMM: fitting GMM (BIC over %d..%d components)...\n', ...
    k_range(1), k_range(end));

gmm_models = cell(1, numel(k_range));
for ki = 1:numel(k_range)
    k = k_range(ki);
    try
        gm = fitgmdist(X, k, ...
            'CovarianceType', o.CovType, ...
            'Replicates',     o.Replicates, ...
            'Options',        statset('MaxIter', 500, 'TolFun', 1e-6), ...
            'RegularizationValue', 1e-6);
        bic_vals(ki)   = gm.BIC;
        gmm_models{ki} = gm;
    catch ME
        warning('classifyTracksGMM: GMM k=%d failed: %s', k, ME.message);
    end
end

[~, best_ki]  = min(bic_vals);
best_k        = k_range(best_ki);
best_gmm      = gmm_models{best_ki};
fprintf('classifyTracksGMM: best k=%d (BIC=%.1f)\n', best_k, bic_vals(best_ki));

labels = cluster(best_gmm, X);   % Nx1, 1..best_k

% ── 5. Sort components by mean α (ascending = slow → fast) ───────────────
mu_alpha = best_gmm.mu(:, 2);   % mean α per component
[~, sort_ord] = sort(mu_alpha, 'ascend');
remap = zeros(1, best_k);
remap(sort_ord) = 1:best_k;
labels_sorted = remap(labels)';

classified_tracks = cell(best_k, 1);
for s = 1:best_k
    classified_tracks{s} = common_track_idx(labels_sorted == s);
end

% ── 6. Optional plot ──────────────────────────────────────────────────────
if o.Plot
    colors = lines(best_k);
    figure;
    hold on;
    for s = 1:best_k
        sel = labels_sorted == s;
        scatter(log10Pc(sel), alpha(sel), 18, colors(s, :), 'filled', ...
            'MarkerFaceAlpha', 0.4, 'DisplayName', sprintf('State %d', s));
    end
    % Draw GMM ellipses (95% confidence)
    theta = linspace(0, 2*pi, 200);
    unit  = [cos(theta); sin(theta)];
    for s = 1:best_k
        orig_comp = find(sort_ord == s);  % original GMM component index
        mu_s  = best_gmm.mu(orig_comp, :);
        sig_s = squeeze(best_gmm.Sigma(:, :, orig_comp));
        if isvector(sig_s)
            sig_s = diag(sig_s);
        end
        % 95% ellipse: chi2inv(0.95,2) ≈ 5.991
        L = chol(sig_s * 5.991, 'lower');
        ellipse = L * unit;
        plot(mu_s(1) + ellipse(1, :), mu_s(2) + ellipse(2, :), '-', ...
            'Color', colors(s, :), 'LineWidth', 1.5);
        plot(mu_s(1), mu_s(2), '+', 'Color', colors(s, :), ...
            'MarkerSize', 12, 'LineWidth', 2);
    end
    xlabel('log_{10}(P_c)');
    ylabel('\alpha (fBM)');
    legend('Location', 'best');
    if ~isempty(o.Title)
        title(o.Title);
    end
    box off;
    set(gca, 'FontSize', 14);

    % BIC plot
    figure;
    plot(k_range, bic_vals, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    xlabel('Number of components');
    ylabel('BIC');
    title('GMM model selection (BIC)');
    set(gca, 'FontSize', 14);
    box off;
end

% ── 7. Pack results ───────────────────────────────────────────────────────
results.track_index        = common_track_idx;
results.log10Pc            = log10Pc;
results.alpha              = alpha;
results.K                  = K;
results.labels             = labels_sorted;
results.gmm                = best_gmm;
results.n_components       = best_k;
results.bic                = bic_vals;
results.classified_tracks  = classified_tracks;

% Print summary
fprintf('\nGMM classification summary (%d tracks, %d states):\n', ...
    numel(common_track_idx), best_k);
fprintf('  %-8s  %-6s  %-12s  %-12s  %-10s\n', ...
    'State', 'N', 'mean log10Pc', 'mean α', 'frac');
for s = 1:best_k
    sel = labels_sorted == s;
    fprintf('  %-8d  %-6d  %-12.2f  %-12.3f  %-10.2f\n', ...
        s, sum(sel), mean(log10Pc(sel)), mean(alpha(sel)), mean(sel));
end
end
