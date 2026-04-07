%% example_smd_to_rl_msd.m
%
% Demonstrates loading a saved SMD object from a .mat file and running a
% Richardson-Lucy decomposition of the MSD distribution via RL_HILO.
%
% CONTEXT
% -------
% A typical upstream workflow uses the SMD class to:
%   1. Load raw images
%   2. Localize single molecules  (smd.localize())
%   3. Link localizations into tracks  (smd.track())
%   4. Optionally compute nuclei masks  (smd.get_roi())
%   5. Cull tracks  (smd.cull_tracks())
%   6. Save the result  (save('my_experiment.mat', 'smd'))
%
% This script picks up from step 6.  It works regardless of whether
% smd.get_roi() was called; see the note in Step 2 below.
%
% PREREQUISITES
% -------------
%   - Run setup_path.m (or add the repo root to the MATLAB path).
%   - The SMD .mat file must be on disk.
%   - Frame rate (Hz) must be known; set SMD_FRAME_RATE below.
%
% OUTPUT
% ------
%   RL_HILO produces two figures per call:
%     Fig 1: van Hove correlation G(r,τ) — data (circles) + RL fit (line)
%     Fig 2: P(MSD) — the Richardson-Lucy deconvolved MSD distribution
%   computed_quantities contains classified tracks, ensemble MSD curves,
%   and power-law fit parameters for each mobility state.

%% 0.  Add pipeline to the path
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setup_path.m'));

%% 1.  Point to your saved SMD .mat file
SMD_MAT_FILE = 'path/to/my_experiment.mat';   % <-- edit this

% Physical parameters — must match the acquisition.
% If smd.frame_rate and smd.pixelsize are already stored in the SMD object
% these values are read automatically in TrajectoryWrapper.fromSMD().
% Only override here if the SMD object is missing those fields.
SMD_FRAME_RATE = 50;    % Hz  (e.g. 50 Hz → 20 ms frame interval)
SMD_PIXEL_SIZE = 0.160; % µm/px  (e.g. 160 nm/px for a typical HILO setup)

%% 2.  Load the SMD object
tmp = load(SMD_MAT_FILE);
% The variable inside the .mat is usually called 'smd'; adjust if different.
smd = tmp.smd;

% --- Handle missing nuclei mask -------------------------------------------
% If smd.get_roi() was never called, smd.ROIs will be empty.
% TrajectoryWrapper.fromSMD() accepts this gracefully: all culled tracks
% (smd.tracks) are imported as-is, with the ROI_ID column taken directly
% from col 10 of the SMD track matrices.
%
% If ROI_IDs in the SMD tracks are all zero (i.e. no ROI assignment),
% getCulledTracks() will return an empty set because it filters col5 > 0.
% In that case, patch the tracks below to assign a dummy ROI_ID of 1:
%
%   for k = 1:numel(smd.tracks)
%       smd.tracks{k}(:, 10) = 1;   % force all tracks into ROI 1
%   end
%
% Comment that block in if you skipped smd.get_roi().
% --------------------------------------------------------------------------

% Patch smd fields that may not exist in older SMD versions
if ~isfield(smd, 'frame_rate') || isempty(smd.frame_rate)
    smd.frame_rate = SMD_FRAME_RATE;
end
if ~isfield(smd, 'pixelsize') || isempty(smd.pixelsize)
    smd.pixelsize = SMD_PIXEL_SIZE;
end
if ~isfield(smd, 'relative_tracks')
    smd.relative_tracks = {};
end
if ~isfield(smd, 'ROIs')
    smd.ROIs = {};
end
if ~isfield(smd, 'fname')
    smd.fname = SMD_MAT_FILE;
end

%% 3.  Wrap the SMD object in a TrajectoryCollection
%
% TrajectoryWrapper.fromSMD() converts the SMD Nx10 track format to the
% pipeline's Nx5 format: [x(µm), y(µm), frame, SNR, ROI_ID].
% FrameInterval is set to 1/smd.frame_rate automatically.

tc = TrajectoryCollection();
tc.addFromSMD(smd, 'FileID', '1', 'Condition', 'Control');

% If you have multiple SMD files (different cells / FOVs), loop here:
%   for k = 1:numel(smd_files)
%       tmp  = load(smd_files{k});
%       tc.addFromSMD(tmp.smd, 'FileID', num2str(k), 'Condition', 'Control');
%   end

tc.summary();

%% 4.  Extract tracks for RL_HILO
%
% RL_HILO (and its sub-functions calc_vanHove, msd_calc, classify_tracks)
% expect tracks in the format:  tracks{1} = cell array of Nx≥2 matrices
%                                with x in col 1, y in col 2, both in µm.
%
% TrajectoryWrapper stores culled tracks in exactly this format after
% fromSMD().  We pull them from the first (and here only) wrapper.

tw     = tc.Wrappers{1};
culled = tw.getCulledTracks();   % cell array of Nx5 matrices [x,y,frame,SNR,ROI_ID]

fprintf('Tracks available for RL analysis: %d\n', numel(culled));

% Wrap in the extra cell layer that RL_HILO / calc_vanHove expect.
tracks_for_rl = {culled};   % tracks_for_rl{1}{k} = k-th track matrix

%% 5.  Richardson-Lucy MSD analysis  (RL_HILO)
%
% RL_HILO deconvolves the van Hove correlation G(r,τ) at a chosen lag time
% to recover the distribution of single-molecule MSDs P(MSD).  It then
% classifies tracks into mobility states and fits ensemble MSD curves to a
% power-law model.
%
% LAGTIME = 1 : use a single-frame displacement (τ = 1 × frame_interval).
%           Increase to e.g. 2 or 3 for slower molecules or noisier data.
%
% TITLE_STR : free-form label used only for figure titles.

LAGTIME   = 1;
TITLE_STR = 'Control';   % e.g. protein name or condition

[M, bins, computed_quantities] = RL_HILO(tracks_for_rl, TITLE_STR, LAGTIME);

%% 6.  Inspect results
%
% computed_quantities fields:
%   .vanHove          : measured G(r,τ) values
%   .Gs               : RL-reconstructed G(r,τ) fit
%   .P1norm           : P(MSD) distribution (normalised)
%   .classified_tracks: {nStates × 1} cell — track indices per mobility state
%   .msd              : [max_lag × nStates] ensemble MSD curves (µm²)
%   .msderr           : corresponding 99 % confidence intervals
%   .fits             : power-law fit parameters per state
%   .fitci            : confidence intervals for fits

n_states = numel(computed_quantities.classified_tracks);
fprintf('\nRL analysis complete: %d mobility state(s) identified.\n', n_states);

dt = tw.getFrameInterval();   % seconds per frame

for s = 1:n_states
    n_tracks = numel(computed_quantities.classified_tracks{s});
    fprintf('  State %d: %d tracks', s, n_tracks);
    if ~isempty(computed_quantities.fits) && numel(computed_quantities.fits) >= s ...
            && ~isempty(computed_quantities.fits{s})
        fp = computed_quantities.fits{s};
        % fp(1) = D (µm²/s^α), fp(2) = offset (localisation noise), fp(3) = α
        fprintf(' | D = %.4f µm²/s^α,  α = %.2f', fp(1) / dt^fp(3), fp(3));
    end
    fprintf('\n');
end

%% 7.  (Optional) Plot ensemble MSD for all states
lag_axis = (1:size(computed_quantities.msd, 1))' * dt;   % convert frames → seconds

figure;
colors = lines(n_states);
hold on;
for s = 1:n_states
    errorbar(lag_axis, computed_quantities.msd(:, s), ...
             computed_quantities.msderr(:, s), ...
             'Color', colors(s,:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('State %d (n=%d)', s, ...
                 numel(computed_quantities.classified_tracks{s})));
end
set(gca, 'XScale', 'log', 'YScale', 'log', ...
         'LineWidth', 1, 'FontSize', 14, 'FontWeight', 'bold');
xlabel('\tau (s)');
ylabel('MSD (\mum^2)');
title([TITLE_STR, ' — Ensemble MSD by state']);
legend('show', 'Location', 'northwest');
box off;

%% 8.  (Optional) Save results
%   save('rl_results_Control.mat', 'M', 'bins', 'computed_quantities');
