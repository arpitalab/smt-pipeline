function [dr2_mean, r_centers, counts] = calc_phi_r(relTracks, lag_frames, varargin)
%CALC_PHI_R  Compute <dR^2>_R — mean squared separation change binned by initial R
%
%  <dR^2>_R is the mean of (Dx)^2 and (Dy)^2 treated as independent 1D samples,
%  binned by the initial scalar separation R = sqrt(dx^2 + dy^2).
%  Each time-pair contributes two samples (one per axis).
%
%  [dr2_mean, r_centers, counts] = calc_phi_r(relTracks, lag_frames)
%  [dr2_mean, r_centers, counts] = calc_phi_r(relTracks, lag_frames, 'NumBins', 40)
%
%  Inputs:
%    relTracks  - cell array of [dx, dy, frame, ROI] matrices (dx,dy in µm)
%    lag_frames - lag in frames (integer)
%
%  Optional name-value:
%    'NumBins'  - number of linear bins (default 40)
%    'MaxR'     - upper bin edge in µm (default 1.0)
%
%  Outputs:
%    dr2_mean  - mean 1D component variance <dR^2>_R per bin (µm^2)
%    r_centers - bin centers (µm)
%    counts    - number of 1D samples per bin

p = inputParser;
addRequired(p, 'relTracks');
addRequired(p, 'lag_frames');
addParameter(p, 'NumBins', 40);
addParameter(p, 'MaxR',    1.0);
parse(p, relTracks, lag_frames, varargin{:});

n_bins = p.Results.NumBins;
max_r  = p.Results.MaxR;
lag    = round(p.Results.lag_frames);

edges = linspace(0, max_r, n_bins + 1);

phi_sum   = zeros(n_bins, 1);
phi_count = zeros(n_bins, 1);

for k = 1:numel(relTracks)
    trk = relTracks{k};
    if size(trk, 1) < 2
        continue;
    end

    frames = trk(:, 3);
    dx     = trk(:, 1);
    dy     = trk(:, 2);
    r      = sqrt(dx.^2 + dy.^2);   % scalar separation, used for binning only

    for i = 1:numel(frames)
        target_frame = frames(i) + lag;
        j = find(frames == target_frame, 1);
        if isempty(j)
            continue;
        end

        r_t  = r(i);
        bin  = find(r_t >= edges(1:end-1) & r_t < edges(2:end), 1);
        if isempty(bin)
            continue;
        end

        % Accumulate x and y components separately (2 samples per time-pair)
        phi_sum(bin)   = phi_sum(bin)   + (dx(j)-dx(i))^2 + (dy(j)-dy(i))^2;
        phi_count(bin) = phi_count(bin) + 2;
    end
end

% Mean 1D component variance per bin
dr2_mean = phi_sum ./ max(phi_count, 1);
dr2_mean(phi_count == 0) = NaN;

r_centers = ((edges(1:end-1) + edges(2:end)) / 2)';
counts    = phi_count;

end
