function results = computePackingFraction(tracks, varargin)
%COMPUTEPACKINGFRACTION  Per-track packing fraction (path tortuosity vs spatial coverage).
%
%   results = computePackingFraction(tracks)
%   results = computePackingFraction(tracks, 'WindowSize', 80, 'MinTrackLength', 100)
%
%   For each track computes:
%
%     pc_window = sum(diff(x).^2 + diff(y).^2) / polyarea(convex_hull)^2
%
%   over sliding windows of WindowSize frames, then takes the median across
%   windows as the per-track aggregate value.
%
%   Interpretation (Renner et al. Biophys J 2017):
%     Pc scales inversely with the confinement area.  Three regimes:
%
%     High pc  —  (1) CONFINED motion: steps and hull area both small,
%                     but hull area² shrinks faster → high Pc.
%                 (2) DIRECTED / sliding motion: steps large, hull area
%                     nearly collinear (→ 0) → Pc highest of all.
%     Low  pc  —  FREE diffusion: hull grows unboundedly → large S² → low Pc.
%
%     To distinguish confinement from directed motion, combine Pc with
%     the hull aspect ratio or straightness index.
%
%   Inputs:
%     tracks  - cell array of N x 2 (or N x 3) trajectory matrices (µm).
%
%   Name-value options:
%     'WindowSize'     - sliding window length in frames (default: 80).
%     'MinTrackLength' - skip tracks shorter than this (default: WindowSize+1).
%     'MinArea'        - skip windows whose convex hull area is below this (µm²)
%                        to avoid division by near-zero for collinear points
%                        (default: 1e-4).
%     'Aggregate'      - 'median' (default) | 'mean' | 'all'
%                        'all' also stores per-window values in results.pc_all.
%     'dt'             - frame interval (s); if provided, also returns
%                        pc_norm = pc / dt (default: []).
%
%   Output struct fields (one entry per qualifying track):
%     .pc           Mx1  per-track aggregate packing fraction (µm^{-2})
%     .pc_norm      Mx1  pc / dt  (if dt provided, else NaN)
%     .n_windows    Mx1  number of valid windows used per track
%     .track_index  Mx1  index into input tracks cell array
%     .pc_all       Mx1  cell of per-window values (only when Aggregate='all')

p = inputParser;
addRequired(p,  'tracks',                    @iscell);
addParameter(p, 'WindowSize',      80,       @isnumeric);
addParameter(p, 'MinTrackLength',  [],       @isnumeric);
addParameter(p, 'MinArea',         1e-4,     @isnumeric);
addParameter(p, 'Aggregate',       'median', @(x) ismember(lower(char(x)), {'median','mean','all'}));
addParameter(p, 'dt',              [],       @isnumeric);
parse(p, tracks, varargin{:});
o = p.Results;

W       = round(o.WindowSize);
Lmin    = o.MinTrackLength;
if isempty(Lmin), Lmin = W + 1; end
agg     = lower(char(o.Aggregate));
want_all = strcmp(agg, 'all');

% Pre-allocate (trim to n_used at the end)
M         = numel(tracks);
pc_out    = nan(M, 1);
pcn_out   = nan(M, 1);
nwin_out  = zeros(M, 1);
idx_out   = zeros(M, 1);
if want_all
    pc_all_cell = cell(M, 1);
end

n_used = 0;

for k = 1:numel(tracks)
    tr = tracks{k};
    L  = size(tr, 1);
    if L < Lmin
        continue;
    end

    x = tr(:, 1);
    y = tr(:, 2);

    win_pc = nan(L - W, 1);
    for m = 1 : L - W
        xw = x(m : m + W - 1);
        yw = y(m : m + W - 1);

        try
            k_hull = convhull(xw, yw);
            A      = polyarea(xw(k_hull), yw(k_hull));
        catch
            continue;   % convhull fails for collinear points
        end

        if A < o.MinArea
            continue;   % near-degenerate (stationary or collinear window)
        end

        ssd        = sum(diff(xw).^2 + diff(yw).^2);
        win_pc(m)  = ssd / A^2;
    end

    valid_wins = win_pc(~isnan(win_pc));
    if isempty(valid_wins)
        continue;
    end

    n_used = n_used + 1;
    idx_out(n_used)  = k;
    nwin_out(n_used) = numel(valid_wins);

    switch agg
        case 'mean'
            pc_out(n_used) = mean(valid_wins);
        otherwise   % 'median' or 'all'
            pc_out(n_used) = median(valid_wins);
    end

    if want_all
        pc_all_cell{n_used} = valid_wins;
    end

    if ~isempty(o.dt)
        pcn_out(n_used) = pc_out(n_used) / o.dt;
    end
end

% Trim to qualifying tracks only
keep = 1:n_used;
results.pc          = pc_out(keep);
results.pc_norm     = pcn_out(keep);
results.n_windows   = nwin_out(keep);
results.track_index = idx_out(keep);
if want_all
    results.pc_all = pc_all_cell(keep);
end

fprintf('computePackingFraction: %d / %d tracks qualified (WindowSize=%d, MinTrackLength=%d)\n', ...
    n_used, numel(tracks), W, Lmin);
end
