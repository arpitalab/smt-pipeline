%% example_msd_and_rl.m
%
% Demonstrates both analysis methods on the same TrajectoryCollection:
%
%   tc.getMSD()              – ensemble-averaged MSD with bootstrap power-law fits
%   tc.getRLDecomposition()  – Richardson-Lucy decomposition of the MSD distribution
%
% getMSD computes the time-averaged mean squared displacement for each track,
% pools them into an ensemble mean, and fits the result to a fractional
% Brownian motion model via bootstrapped lsqnonlin.
%
% getRLDecomposition deconvolves the van Hove displacement distribution at a
% chosen lag time to recover P(MSD) — the distribution of single-molecule MSDs
% — and classifies tracks into mobility states.
%
% PREREQUISITES
% -------------
%   - Run setup_path.m (or add the repo root to the MATLAB path).
%   - CSV files with columns: x, y, frame, track_id (and optionally snr).
%   - A TOML parameter file specifying at minimum pixel_size, frame_interval,
%     and exposure_time_s (needed by getMSD to compute the motion-blur
%     correction factor frac = exposure_time / frame_interval).

%% 0.  Path
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setup_path.m'));

%% 1.  Build a TrajectoryCollection
%
% Edit DATA_DIR and TOML_FILE for your experiment.

DATA_DIR  = 'path/to/csv/files';
TOML_FILE = 'path/to/experiment.toml';
% Minimal experiment.toml:
%   [acquisition]
%   pixel_size      = 0.160   # µm/px
%   frame_interval  = 0.020   # s  (50 Hz)
%   exposure_time_s = 0.010   # s  (10 ms exposure)

tc = TrajectoryCollection();
tc.ReadParams(TOML_FILE);

files = dir(fullfile(DATA_DIR, '*.csv'));
for k = 1:numel(files)
    tc.addFromFile(fullfile(files(k).folder, files(k).name), ...
                   'FileID',    num2str(k), ...
                   'Condition', 'Control');
end
tc.summary();

% =========================================================================
%% PART A: Ensemble MSD  (tc.getMSD)
% =========================================================================

%% 2.  Run ensemble MSD
%
% TrackType  : 'culled'   – tracks inside nucleus mask (default)
%              'relative' – after global cell-motion removal
%              'raw'      – unfiltered tracks
% MaxLag     : maximum lag in frames over which MSD is computed and fit
% NumBootstrap: number of bootstrap resamples for confidence intervals
% ExposureTime: camera exposure in seconds; used to compute the motion-blur
%               correction frac = ExposureTime / FrameInterval.
%               Read automatically from Params.exposure_time_s if set in
%               the TOML; pass explicitly to override.

msd_results = tc.getMSD( ...
    'TrackType',    'culled', ...
    'MaxLag',       25,       ...
    'NumBootstrap', 50,       ...
    'ExposureTime', 0.01);    % override if not in TOML

%% 3.  Inspect MSD results
%
% msd_results fields:
%   .ensemble_mean  – [MaxLag × 1] ensemble-averaged MSD (µm²)
%   .lag_axis       – [MaxLag × 1] lag times in seconds
%   .lag_frames     – [MaxLag × 1] integer lag indices
%   .boot_means     – [MaxLag × NumBootstrap] per-bootstrap ensemble means
%   .boot_fits      – [NumBootstrap × 3] fit parameters [G, sigma_sq, alpha]
%   .dt             – frame interval (s)
%   .frac           – exposure / frame_interval
%   .n_tracks       – number of tracks used

fprintf('\nEnsemble MSD: %d tracks, dt = %.4f s, frac = %.3f\n', ...
    msd_results.n_tracks, msd_results.dt, msd_results.frac);

% Median fit parameters across bootstrap samples
valid_fits = msd_results.boot_fits(~any(isnan(msd_results.boot_fits), 2), :);
if ~isempty(valid_fits)
    med_G     = median(valid_fits(:, 1));
    med_sig2  = median(valid_fits(:, 2));
    med_alpha = median(valid_fits(:, 3));
    ci_alpha  = prctile(valid_fits(:, 3), [2.5, 97.5]);
    fprintf('Median fit:  G = %.4f µm²/s^α,  σ² = %.4f µm²,  α = %.2f  [%.2f – %.2f]\n', ...
        med_G, med_sig2, med_alpha, ci_alpha(1), ci_alpha(2));
end

%% 4.  Plot ensemble MSD with bootstrap confidence band
tau   = msd_results.lag_axis;
emsd  = msd_results.ensemble_mean;
bmean = msd_results.boot_means;

ci_lo = prctile(bmean, 2.5,  2);
ci_hi = prctile(bmean, 97.5, 2);

figure('Name', 'Ensemble MSD');
patch([tau; flipud(tau)], [ci_lo; flipud(ci_hi)], [0.7 0.85 1], ...
      'EdgeColor', 'none', 'FaceAlpha', 0.6);
hold on;
plot(tau, emsd, 'b-', 'LineWidth', 2);

% Overlay median fit curve
if ~isempty(valid_fits)
    lag_vec = msd_results.lag_frames;
    msd_fit = my_fun([med_G, med_sig2, med_alpha], lag_vec, emsd, ...
                     msd_results.dt, msd_results.frac);
    % my_fun returns residuals; reconstruct the model curve
    G = med_G; a = med_alpha; s2 = med_sig2;
    dt = msd_results.dt; fr = msd_results.frac;
    b_vec = (abs(1 + fr./lag_vec).^(2+a) + abs(1-fr./lag_vec).^(2+a) - 2) ./ (fr./lag_vec).^2;
    fit_curve = G/((1+a)*(2+a)) * ((dt.*lag_vec).^a .* b_vec - 2*(fr*dt)^a) + 2*s2;
    plot(tau, fit_curve, 'r--', 'LineWidth', 1.5);
    legend({'95% CI (bootstrap)', 'Ensemble mean', 'fBM fit (median)'}, 'Location', 'northwest');
else
    legend({'95% CI (bootstrap)', 'Ensemble mean'}, 'Location', 'northwest');
end

set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 14, 'LineWidth', 1);
xlabel('\tau (s)');
ylabel('MSD (\mum^2)');
title('Ensemble MSD — Control');
box off;

% =========================================================================
%% PART B: Richardson-Lucy MSD decomposition  (tc.getRLDecomposition)
% =========================================================================

%% 5.  Run RL decomposition
%
% LagTime : frame lag for the van Hove correlation.
%           Use a lag where mobile and immobile populations are well separated
%           (typically 1–5 for slow imaging, higher for fast).
% String  : label for figure titles only.
%
% Produces two figures automatically:
%   Fig A: G(r,τ) — van Hove correlation, measured vs. RL fit
%   Fig B: P(MSD) — deconvolved MSD probability distribution

rl_results = tc.getRLDecomposition('LagTime', 4, 'String', 'H2B Control');

%% 6.  Inspect RL results
%
% rl_results fields:
%   .M                   – MSD bin centers for P(MSD)
%   .computed_quantities – struct with full RL output:
%       .vanHove          – measured G(r,τ)
%       .Gs               – RL-reconstructed G(r,τ)
%       .P1norm           – P(MSD) distribution (normalised)
%       .classified_tracks – {nStates × 1} cell of track index lists
%       .msd              – [max_lag × nStates] ensemble MSD per state
%       .msderr           – 99% CI on state MSDs
%       .fits             – power-law fit params per state

cq       = rl_results.computed_quantities;
n_states = numel(cq.classified_tracks);
fprintf('\nRL decomposition: %d mobility state(s) at lag %d.\n', ...
    n_states, cq.lagtime);
for s = 1:n_states
    fprintf('  State %d: %d tracks\n', s, numel(cq.classified_tracks{s}));
end

%% 7.  Plot per-state ensemble MSD from RL classification
%
% The RL method classifies each track into a mobility state.
% Here we plot the ensemble MSD for each state on the same axes.

dt     = tc.Wrappers{1}.getFrameInterval();
n_lags = size(cq.msd, 1);
tau_rl = (1:n_lags)' * dt;

figure('Name', 'RL per-state MSD');
colors = lines(n_states);
hold on;
for s = 1:n_states
    errorbar(tau_rl, cq.msd(:, s), cq.msderr(:, s), ...
             'Color', colors(s, :), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('State %d  (n=%d)', s, numel(cq.classified_tracks{s})));
end
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 14, 'LineWidth', 1);
xlabel('\tau (s)');
ylabel('MSD (\mum^2)');
title('Per-state MSD — H2B Control');
legend('show', 'Location', 'northwest');
box off;

%% 8.  (Optional) Save
%   tc.save('H2B_Control.mat');
