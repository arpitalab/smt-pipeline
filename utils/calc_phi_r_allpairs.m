function [dr2_mean, r_centers, counts] = calc_phi_r_allpairs(source, lag_frames, varargin)
%CALC_PHI_R_ALLPAIRS  Compute <dR^2>_R over all within-nucleus track pairs
%
%  Loops over all nuclei (ROIs), pairing every track in the same nucleus,
%  finds frames where both tracks and their lag counterparts exist, and
%  accumulates squared separation change binned by initial separation R.
%
%  ROI IDs are only treated as unique within a single TrajectoryWrapper.
%  When a TrajectoryCollection is passed, each wrapper is processed
%  separately to avoid ROI ID collisions across files.
%
%  Usage:
%    [dr2_mean, r_centers, counts] = calc_phi_r_allpairs(tc, lag_frames)
%    [dr2_mean, r_centers, counts] = calc_phi_r_allpairs(tw, lag_frames)
%    [dr2_mean, r_centers, counts] = calc_phi_r_allpairs(tc, lag_frames, 'MaxPairs', 200)
%
%  Inputs:
%    source     - TrajectoryCollection or TrajectoryWrapper
%    lag_frames - lag in frames (integer)
%
%  Optional name-value:
%    'NumBins'  - number of linear bins from 0 to MaxR (default 40)
%    'MaxR'     - upper bin edge in µm (default 1.0)
%    'MaxPairs' - max pairs sampled per nucleus; 0 = use all (default 0)
%
%  Outputs:
%    dr2_mean  - mean 1D component variance <dR^2>_R per bin (µm^2)
%    r_centers - bin centers (µm)
%    counts    - number of 1D samples per bin (each frame-pair -> 2 samples)

p = inputParser;
addRequired(p, 'source');
addRequired(p, 'lag_frames');
addParameter(p, 'NumBins',     40);
addParameter(p, 'MaxR',        1.0);
addParameter(p, 'MaxPairs',    0);
addParameter(p, 'MaxWrappers', 0);   % 0 = use all
parse(p, source, lag_frames, varargin{:});

n_bins    = p.Results.NumBins;
max_r     = p.Results.MaxR;
lag       = round(p.Results.lag_frames);
max_pair  = p.Results.MaxPairs;
max_wrap  = p.Results.MaxWrappers;

% Collect wrappers
if isa(source, 'TrajectoryCollection')
    wrappers = source.Wrappers;
elseif isa(source, 'TrajectoryWrapper')
    wrappers = {source};
else
    error('calc_phi_r_allpairs: source must be a TrajectoryWrapper or TrajectoryCollection');
end

bin_w     = max_r / n_bins;
phi_sum   = zeros(n_bins, 1);
phi_count = zeros(n_bins, 1);

n_wrap = numel(wrappers);
if max_wrap > 0
    n_wrap = min(max_wrap, n_wrap);
end

for wi = 1:n_wrap
    tracks = wrappers{wi}.getCulledTracks();
    if isempty(tracks)
        continue;
    end
    [phi_sum, phi_count] = accumulate_wrapper(tracks, lag, n_bins, bin_w, max_pair, phi_sum, phi_count);
end

dr2_mean = phi_sum ./ max(phi_count, 1);
dr2_mean(phi_count == 0) = NaN;

edges     = linspace(0, max_r, n_bins + 1);
r_centers = ((edges(1:end-1) + edges(2:end)) / 2)';
counts    = phi_count;

end

% -------------------------------------------------------------------------
function [phi_sum, phi_count] = accumulate_wrapper(tracks, lag, n_bins, bin_w, max_pair, phi_sum, phi_count)

N = numel(tracks);
trk(N) = struct('frames', [], 'x', [], 'y', [], 'roi', []);
for k = 1:N
    m             = tracks{k};
    [f, ord]      = sort(m(:,3));
    trk(k).frames = f;
    trk(k).x      = m(ord, 1);
    trk(k).y      = m(ord, 2);
    trk(k).roi    = m(1, 5);
end

rois       = [trk.roi];
unique_roi = unique(rois);

for ri = 1:numel(unique_roi)
    idx = find(rois == unique_roi(ri));
    n   = numel(idx);
    if n < 2
        continue;
    end

    % All pairs within this nucleus, with optional subsampling
    [aa, bb] = find(tril(ones(n), -1));
    n_pairs  = numel(aa);
    if max_pair > 0 && n_pairs > max_pair
        sel = randperm(n_pairs, max_pair);
        aa  = aa(sel);
        bb  = bb(sel);
    end

    for pi = 1:numel(aa)
        ia = idx(aa(pi));
        ib = idx(bb(pi));

        fa = trk(ia).frames;
        fb = trk(ib).frames;

        common  = intersect(fa, fb);
        valid_t = common(ismember(common + lag, common));
        if isempty(valid_t)
            continue;
        end

        [~, ia_t]  = ismember(valid_t,       fa);
        [~, ia_tl] = ismember(valid_t + lag, fa);
        [~, ib_t]  = ismember(valid_t,       fb);
        [~, ib_tl] = ismember(valid_t + lag, fb);

        rel_x_t  = trk(ia).x(ia_t)  - trk(ib).x(ib_t);
        rel_y_t  = trk(ia).y(ia_t)  - trk(ib).y(ib_t);
        rel_x_tl = trk(ia).x(ia_tl) - trk(ib).x(ib_tl);
        rel_y_tl = trk(ia).y(ia_tl) - trk(ib).y(ib_tl);

        R_t  = sqrt(rel_x_t.^2 + rel_y_t.^2);
        dRx2 = (rel_x_tl - rel_x_t).^2;
        dRy2 = (rel_y_tl - rel_y_t).^2;

        bin_idx = ceil(R_t / bin_w);
        valid   = bin_idx >= 1 & bin_idx <= n_bins;
        bin_idx = bin_idx(valid);
        if isempty(bin_idx)
            continue;
        end
        dr2_vec = dRx2(valid) + dRy2(valid);

        phi_sum   = phi_sum   + accumarray(bin_idx, dr2_vec, [n_bins, 1]);
        phi_count = phi_count + accumarray(bin_idx, 2,       [n_bins, 1]);
    end
end

end
