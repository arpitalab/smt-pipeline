% setup_path  Add all package folders to the MATLAB path.
%
% Run once per MATLAB session before using the smt-pipeline:
%
%   run('path/to/smt-pipeline/setup_path.m')
%
% After this, TrajectoryWrapper, TrajectoryCollection, TrajectoryAdapter,
% and TrackUtils are available without further path manipulation.

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, 'utils'));
addpath(fullfile(here, 'utils', 'RL_analysis'));
addpath(fullfile(here, 'third_party', 'matlab-toml'));

fprintf('smt-pipeline path configured.\n');
