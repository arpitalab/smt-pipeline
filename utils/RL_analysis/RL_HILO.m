function [M, bins, computed_quantities] = RL_HILO(tracks, title_str, lagtime, varargin)
% RL_HILO  Richardson-Lucy decomposition of the MSD distribution.
%
%   [M, bins, computed_quantities] = RL_HILO(tracks, title_str, lagtime)
%   [M, bins, computed_quantities] = RL_HILO(tracks, title_str, lagtime, ...
%                                       'dt', 0.02, 'ExposureFraction', 0.5)
%
%   Inputs:
%     tracks        - Cell array in Kaustubh format: tracks{1}{itrack} is Nx2 (x,y µm)
%     title_str     - Label for figure titles
%     lagtime       - Frame lag for van Hove / track classification
%
%   Name-value options:
%     'dt'               - Frame interval in seconds (default: 0.02 s).
%                          Used for MSD axis and fBM fitting.
%     'ExposureFraction' - Te/dt ratio for fBM model (default: 1).
%     'MaxLag'           - Maximum lag for MSD calculation (default: 25).
%     'MinGroupSize'     - Minimum tracks in a state to fit fBM (default: 50).
%     'SubtrackLength'   - Subtrack length for population fBM MLE (default: 20).
%     'CIMethod'         - CI method for fitFBM_MLE (default: 'none').
%
%   Outputs:
%     M                   - Vector of MSD values (RL deconvolution result)
%     bins                - Displacement values for van Hove correlation
%     computed_quantities - Struct with full results (see below)
%
%   computed_quantities fields:
%     .ident             - title_str
%     .lagtime           - lagtime used
%     .dt                - frame interval (s)
%     .frac              - ExposureFraction used
%     .x                 - bins (displacements)
%     .vanHove           - empirical van Hove correlation
%     .P1norm            - RL-deconvolved P(MSD)
%     .Gs                - estimated van Hove from RL fit
%     .classified_tracks - cell array of track index vectors, one per state
%     .msd               - MaxLag x nStates ensemble MSD matrix
%     .msderr            - MaxLag x nStates MSD error matrix
%     .fbm               - nStates x 1 struct array: per-state fBM MLE results
%                          fields: .K, .alpha, .Ka, .sigma, .n_tracks, .converged
%                          (empty if state has < MinGroupSize tracks)

p = inputParser;
addParameter(p, 'dt',               0.02,     @isnumeric);
addParameter(p, 'ExposureFraction', 1,        @isnumeric);
addParameter(p, 'MaxLag',           25,       @isnumeric);
addParameter(p, 'MinGroupSize',     50,       @isnumeric);
addParameter(p, 'SubtrackLength',   20,       @isnumeric);
addParameter(p, 'CIMethod',         'none',   @ischar);
parse(p, varargin{:});
o = p.Results;

dt      = o.dt;
frac    = o.ExposureFraction;
max_lag = round(o.MaxLag);

fprintf('RL_HILO: dt=%.4f s, ExposureFraction=%.3f\n', dt, frac);

% --- Van Hove and RL deconvolution ---
[~, vanHove, bins] = calc_vanHove(tracks, lagtime, 1);
[M, Gs, P1norm, ~] = RL_psd_tmp(bins, vanHove, 50000);

figure;
plot(bins, vanHove, 'ko', 'markersize', 8);
hold on;
set(gca, 'yscale', 'log');
plot(bins, Gs, 'r', 'linewidth', 2);
set(gca, 'linewidth', 1, 'fontsize', 24, 'fontweight', 'bold');
xlabel('r [\mum]');
ylabel('G(r,\tau)');
axis([0 1 1e-4 100]);
axis square;
box off;
title(title_str);

figure;
semilogx(M, P1norm, 'k', 'linewidth', 2);
title(title_str);
axis([1e-3 2e-1 0 200]);
xlabel('MSD [\mum^2]');
ylabel('P(MSD)');
set(gca, 'linewidth', 1, 'fontsize', 24, 'fontweight', 'bold');
box off;

% --- Track classification ---
classified_tracks = classify_tracks(tracks, P1norm, M, lagtime, 1, title_str);
n_states = numel(classified_tracks);

% --- Ensemble MSD per state ---
[x_lag, msd, msd_err] = msd_calc(classified_tracks, tracks, max_lag);

% --- Per-state fBM MLE (replaces lsqnonlin + my_fun) ---
fbm_results(n_states) = struct('K', NaN, 'alpha', NaN, 'Ka', NaN, ...
    'sigma', NaN, 'n_tracks', 0, 'converged', false, 'CI', struct());

all_tracks_flat = tracks{1};  % unwrap Kaustubh format

for s = 1:n_states
    idx = classified_tracks{s};
    n_s = numel(idx);

    if n_s < o.MinGroupSize
        fprintf('RL_HILO: state %d — %d tracks (< MinGroupSize=%d), skipping fBM fit.\n', ...
            s, n_s, o.MinGroupSize);
        fbm_results(s).n_tracks = n_s;
        continue;
    end

    state_tracks = all_tracks_flat(idx);

    fprintf('RL_HILO: fitting fBM MLE for state %d (%d tracks)...\n', s, n_s);
    try
        [fp, ~, ci] = fitFBM_MLE(state_tracks, dt, ...
            'ExposureFraction', frac, ...
            'SubtrackLength',   o.SubtrackLength, ...
            'CIMethod',         o.CIMethod, ...
            'Verbose',          false);
        fbm_results(s).K         = fp.K;
        fbm_results(s).alpha     = fp.alpha;
        fbm_results(s).Ka        = fp.Ka;
        fbm_results(s).sigma     = fp.sigma;
        fbm_results(s).n_tracks  = n_s;
        fbm_results(s).converged = true;
        fbm_results(s).CI        = ci;
        fprintf('  state %d: alpha=%.3f, K=%.4g µm²/s^alpha, sigma=%.4f µm\n', ...
            s, fp.alpha, fp.K, fp.sigma);
    catch ME
        warning('RL_HILO: fBM fit failed for state %d: %s', s, ME.message);
        fbm_results(s).n_tracks = n_s;
    end
end

% --- Pack output ---
computed_quantities.ident             = title_str;
computed_quantities.lagtime           = lagtime;
computed_quantities.dt                = dt;
computed_quantities.frac              = frac;
computed_quantities.x                 = bins;
computed_quantities.vanHove           = vanHove;
computed_quantities.P1norm            = P1norm;
computed_quantities.Gs                = Gs;
computed_quantities.classified_tracks = classified_tracks;
computed_quantities.msd               = msd;
computed_quantities.msderr            = msd_err;
computed_quantities.msd_lags          = x_lag;
computed_quantities.fbm               = fbm_results;

end
