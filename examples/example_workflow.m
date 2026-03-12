%% example_workflow.m
%
% Minimal end-to-end demonstration of the smt-pipeline.
%
% Prerequisites:
%   1. Run setup_path.m (or add the repo root to your MATLAB path).
%   2. Install pEMv2 and add it to the path (see DEPENDENCIES.md).
%   3. Replace DATA_DIR and TOML_FILE with your own paths.
%
% The pipeline expects CSV files with at least these columns:
%   x, y      : localisation coordinates (pixels)
%   frame     : frame index (1-based integer)
%   track_id  : integer track identifier
%   snr       : signal-to-noise ratio (optional but recommended)

%% 0. Configure paths
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setup_path.m'));

DATA_DIR  = 'path/to/csv/files';   % folder containing one CSV per cell
TOML_FILE = 'path/to/experiment.toml';  % parameter file (see below)

%% 1. Create a TrajectoryCollection and load parameters
%
% A minimal experiment.toml looks like:
%
%   [acquisition]
%   pixel_size     = 0.165   # µm per pixel
%   frame_interval = 0.02    # seconds between frames
%
%   [culling]
%   min_track_length = 5
%
tc = TrajectoryCollection();
tc.ReadParams(TOML_FILE);

%% 2. Add CSV files
%
% Each addFromFile call:
%   - Loads the CSV, scales x/y to µm using PixelSize
%   - Culls short/low-quality tracks
%   - Removes global cell motion via maximum spanning forest
%   - Stores the resulting TrajectoryWrapper internally

files = dir(fullfile(DATA_DIR, '*.csv'));
if isempty(files)
    error('No CSV files found in %s', DATA_DIR);
end

for k = 1:numel(files)
    tc.addFromFile(fullfile(files(k).folder, files(k).name), ...
                   'FileID',    num2str(k), ...
                   'Condition', 'Control');
    fprintf('Loaded %d / %d: %s\n', k, numel(files), files(k).name);
end

%% 3. Summary statistics
%
% Prints per-file track counts, mean track length, and fraction of
% mobile vs. immobile particles.

tc.summary();

%% 4. MSD analysis
%
% Computes ensemble-averaged MSD for each condition.
% Results are stored in tc.MSDResults and can be plotted with tc.plotMSD().

tc.computeMSD();

%% 5. Bayesian diffusivity classification (pEMv2)
%
% Requires pEMv2 on the MATLAB path (see DEPENDENCIES.md).
% Classifies each track as mobile/immobile and estimates diffusion
% coefficients using a permutation expectation-maximisation approach.

tc.getBayesianDiffusivity('Condition', 'Control');

%% 6. Lifetime analysis
%
% Computes bound-particle dwell-time distributions.

tc.computeLifetimes('Condition', 'Control');

%% Done
fprintf('\nWorkflow complete.\n');
