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

valid_fits = msd_results.boot_fits(~any(isnan(msd_results.boot_fits), 2), :);
if ~isempty(valid_fits)
    med_alpha = median(valid_fits(:, 3));
    ci_alpha  = prctile(valid_fits(:, 3), [2.5, 97.5]);
    fprintf('Median fit:  G = %.4f µm²/s^α,  σ² = %.4f µm²,  α = %.2f  [%.2f–%.2f]\n', ...
        median(valid_fits(:,1)), median(valid_fits(:,2)), med_alpha, ci_alpha(1), ci_alpha(2));
end

%% 4.  Plot ensemble MSD
%
% tc.plotMSD() produces one figure:
%   - shaded 95% bootstrap CI band
%   - ensemble mean line
%   - median fBM fit curve with α in the legend
%
% Optional: 'ShowBootstraps', true overlays individual bootstrap fit curves.

tc.plotMSD('Title', 'H2B Control', 'Color', [0 0.4 0.8]);
% tc.plotMSD('Title', 'H2B Control', 'ShowBootstraps', true);

% =========================================================================
%% PART B: Richardson-Lucy MSD decomposition  (tc.getRLDecomposition)
% =========================================================================

%% 5.  Run RL decomposition
%
% LagTime : frame lag for the van Hove correlation.
%           Use a lag where mobile and immobile populations are well separated
%           (typically 1–5 for slow imaging, higher for fast).
% String  : label for figure titles only.

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

cq = rl_results.computed_quantities;
fprintf('\nRL decomposition: %d mobility state(s) at lag %d.\n', ...
    numel(cq.classified_tracks), cq.lagtime);
for s = 1:numel(cq.classified_tracks)
    fprintf('  State %d: %d tracks\n', s, numel(cq.classified_tracks{s}));
end

%% 7.  Plot RL results
%
% tc.plotRLDecomposition() produces three figures:
%   Fig 1 – van Hove G(r,τ): measured data vs. RL fit
%   Fig 2 – P(MSD): Richardson-Lucy deconvolved MSD distribution
%   Fig 3 – Per-state ensemble MSD with 99% CI error bars

tc.plotRLDecomposition('Title', 'H2B Control');

%% 8.  (Optional) Save
%   tc.save('H2B_Control.mat');
