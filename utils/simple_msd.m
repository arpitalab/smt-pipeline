function [msd, num_lags_tracks, lags] = simple_msd(pos, t)
% SIMPLE_MSD  Per-track MSD calculation supporting gapped (NaN-filled) trajectories.
%
%   [msd, num_lags_tracks, lags] = simple_msd(pos, t)
%
%   Inputs:
%     pos - Nx2 position matrix [x, y] (µm); NaN entries are handled gracefully
%     t   - Nx1 time vector (uniform or non-uniform, in frame units or seconds)
%
%   Outputs:
%     msd             - sum of squared displacements at each lag (not yet divided
%                       by num_lags_tracks; divide outside to get ensemble mean)
%     num_lags_tracks - number of valid displacement pairs at each lag
%     lags            - lag values in the same units as t

delta_t = diff(t);
x = pos(:, 1);
y = pos(:, 2);

if min(delta_t) == max(delta_t)
    % Uniform time — standard lag loop
    max_lag_track = length(t) - 1;
    tmp_msd           = zeros(max_lag_track, 1);
    tmp_num_lags      = zeros(max_lag_track, 1);
    for lagtime = 1:max_lag_track
        r = (x(lagtime+1:end) - x(1:end-lagtime)).^2 + ...
            (y(lagtime+1:end) - y(1:end-lagtime)).^2;
        tmp_msd(lagtime)      = sum(r, 'omitnan');
        tmp_num_lags(lagtime) = sum(~isnan(r));
    end
    msd             = tmp_msd;
    num_lags_tracks = tmp_num_lags;
    lags            = (delta_t(1):delta_t(1):max_lag_track * delta_t(1))';
else
    % Non-uniform time — pairwise distance approach
    interval = min(delta_t);
    DeltaT   = t(end) - t(1);
    lag_idx  = round((t(1:2:end) - t(1)) / (2 * interval));
    dt_pairs = pdist(lag_idx);
    maxlag   = round(DeltaT / interval) - 1;
    Dlag     = min(4000, round(maxlag));
    D        = pdist(pos(:, 1:2), 'squaredeuclidean');
    msd             = zeros(Dlag, 1);
    num_lags_tracks = zeros(Dlag, 1);
    for jj = 1:Dlag
        ix               = (dt_pairs == jj);
        msd(jj)          = sum(D(ix));
        num_lags_tracks(jj) = sum(ix);
    end
    lags = (1:Dlag)' * interval;
end
end
