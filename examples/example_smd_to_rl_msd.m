%% example_smd_to_rl_msd.m
%
% Load a saved SMD object and run a Richardson-Lucy MSD decomposition.
%
% UPSTREAM WORKFLOW
% -----------------
% The SMD object is produced by:
%   smd.localize()      % detect single molecules in raw images
%   smd.track()         % link detections into trajectories
%   smd.get_roi()       % compute nuclei masks  (optional)
%   smd.cull_tracks()   % filter short / out-of-nucleus tracks
%   save('my_experiment.mat', 'smd')
%
% This script picks up from the saved .mat file.
%
% PREREQUISITES
% -------------
%   - Run setup_path.m (or add the repo root to the MATLAB path).
%   - Frame rate (Hz) must be stored in smd.frame_rate; set the fallback
%     value SMD_FRAME_RATE below if the field is missing.
%
% OUTPUT
% ------
%   tc.getRLDecomposition() calls RL_HILO internally and produces two figures:
%     Fig 1: van Hove correlation G(r,τ) — data + RL fit
%     Fig 2: P(MSD) — Richardson-Lucy deconvolved MSD distribution
%   Results are stored in tc.RLResults.

%% 0.  Add pipeline to the path
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setup_path.m'));

%% 1.  Load the SMD object
SMD_MAT_FILE   = 'path/to/my_experiment.mat';   % <-- edit this
SMD_FRAME_RATE = 50;    % Hz — fallback if smd.frame_rate is missing

tmp = load(SMD_MAT_FILE);
smd = tmp.smd;   % adjust variable name if different inside the .mat

% Patch fields that may be absent in older SMD objects
if ~isfield(smd, 'frame_rate') || isempty(smd.frame_rate)
    smd.frame_rate = SMD_FRAME_RATE;
end
if ~isfield(smd, 'relative_tracks'), smd.relative_tracks = {}; end
if ~isfield(smd, 'ROIs'),            smd.ROIs = {};            end
if ~isfield(smd, 'fname'),           smd.fname = SMD_MAT_FILE; end

% If smd.get_roi() was never called, ROI_IDs in the track matrices (col 10)
% may be zero, which would cause getCulledTracks() to return nothing.
% Uncomment the block below to assign a dummy ROI_ID and proceed anyway:
%
%   for k = 1:numel(smd.tracks)
%       smd.tracks{k}(:, 10) = 1;
%   end

%% 2.  Build a TrajectoryCollection from the SMD object
%
% fromSMD() converts the SMD Nx10 track format to the pipeline's Nx5
% format [x(µm), y(µm), frame, SNR, ROI_ID] and sets FrameInterval
% automatically from smd.frame_rate.

tc = TrajectoryCollection();
tc.addFromSMD(smd, 'FileID', '1', 'Condition', 'Control');

% Multiple cells / FOVs — just loop:
%   smd_files = {'cell1.mat', 'cell2.mat', ...};
%   for k = 1:numel(smd_files)
%       tmp = load(smd_files{k});
%       tc.addFromSMD(tmp.smd, 'FileID', num2str(k), 'Condition', 'Control');
%   end

tc.summary();

%% 3.  Richardson-Lucy MSD decomposition
%
% LagTime : frame lag at which the van Hove correlation is computed.
%           Lag 4 (= 4 × frame_interval seconds) is a common starting
%           point; adjust based on your molecule's timescale.
% String  : label used in figure titles only.

results = tc.getRLDecomposition('LagTime', 4, 'String', 'H2B');

%% 4.  (Optional) Save results
%   tc.save('H2B_Control.mat');
