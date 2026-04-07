%% example_smtdatabase.m
%
% Demonstrates using SMTDatabase to register experiments, log analyses,
% browse results, and reload saved collections across sessions.
%
% SMTDatabase is a lightweight SQLite registry (no Database Toolbox needed,
% requires MATLAB R2022b+).  It stores experiment metadata, links FOVs to
% condition groups ("collections"), and records analysis results so you can
% track what has been run and quickly reload saved .mat files.
%
% PREREQUISITES
% -------------
%   - Run setup_path.m (or add the repo root to the MATLAB path).
%   - Have at least one TrajectoryCollection built from CSV or SMD files.
%     See example_workflow.m or example_smd_to_rl_msd.m.

%% 0.  Add pipeline to the path
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setup_path.m'));

% =========================================================================
%% PART 0 (OPTIONAL): UPSTREAM SMD PROCESSING
%
% Skip this section if you already have CSV trajectory files or a saved
% SMD .mat file — jump straight to Part A.
%
% This section shows how to run the full SMD localization-and-tracking
% pipeline on a folder of TIFF or ND2 movies, save each result, and
% collect the processed SMD objects into a TrajectoryCollection.
%
% The SMD class must be on the MATLAB path (it lives in sr_tracking/@SMD/).
% =========================================================================

%% 0a.  Set acquisition parameters (shared across all files in the batch)

IMAGE_DIR    = 'path/to/tiff_or_nd2_files';   % <-- folder containing raw movies
EXPOSURE_S   = 0.010;    % camera exposure time in seconds (e.g. 10 ms)
FRAME_RATE   = 5;        % acquisition frame rate in Hz  (e.g. 5 Hz = 200 ms interval)
SAVE_DIR     = 'path/to/smd_output';          % where to write processed .mat files
CONDITION    = 'Control';

% Create output folder if it does not exist
if ~isfolder(SAVE_DIR), mkdir(SAVE_DIR); end

%% 0b.  Process each movie file

img_files = [dir(fullfile(IMAGE_DIR, '*.tif')); ...
             dir(fullfile(IMAGE_DIR, '*.nd2'))];

tc_from_smd = TrajectoryCollection();

for k = 1:numel(img_files)
    fname   = img_files(k).name;
    fdir    = img_files(k).folder;
    fprintf('\n[%d/%d] Processing %s\n', k, numel(img_files), fname);

    % --- Construct SMD object -------------------------------------------
    smd = SMD(fdir, fname, EXPOSURE_S, FRAME_RATE);

    % Detection / localization settings
    smd.detection_method        = 'wavelet'; % 'wavelet' (default) or 'llr'
    smd.localization_threshold  = 1.75;      % wavelet threshold; lower = more spots
    smd.localization_box        = 3;         % half-width of fitting box (pixels)
    smd.pixelsize               = 0.160;     % µm/px

    % For LLR detection instead of wavelet, uncomment and tune:
    %   smd.detection_method = 'llr';
    %   smd.noise_model      = 'gaussian';   % 'gaussian' or 'poisson'
    %   smd.psf_sigma        = 1.3;          % PSF std dev in pixels
    %   smd.llr_params.t     = 20.0;         % score threshold

    % Tracking parameters
    smd.tracking_params.trackingRadius     = 2;   % max displacement per frame (pixels)
    smd.tracking_params.gapFrames          = 3;   % max gap (missing frames) to bridge
    smd.tracking_params.minLengthBeforeGap = 4;   % min run length before first gap

    % --- Run pipeline ---------------------------------------------------
    smd.localize();   % detect spots and fit sub-pixel positions

    % Nucleus / ROI mask (optional but recommended for nuclear proteins).
    % Requires a matching fluorescence image or BF image in the same folder.
    % Comment out if no mask is available; addFromSMD handles empty ROIs.
    try
        smd.get_roi();
    catch ME
        warning('get_roi failed for %s: %s', fname, ME.message);
    end

    smd.track();        % link localizations into trajectories
    smd.cull_tracks();  % remove short / out-of-nucleus tracks

    % --- Save the SMD object -------------------------------------------
    [~, stem] = fileparts(fname);
    save_path = fullfile(SAVE_DIR, [stem '_smd.mat']);
    save(save_path, 'smd');
    fprintf('    Saved → %s\n', save_path);

    % --- Add to TrajectoryCollection -----------------------------------
    tc_from_smd.addFromSMD(smd, 'FileID', num2str(k), 'Condition', CONDITION);
end

tc_from_smd.summary();

% tc_from_smd is now equivalent to a collection built from CSV files.
% Continue with Part A using tc_from_smd in place of tc, or proceed
% directly to analysis and database registration below.

%% 1.  Open (or create) the database
%
% If the .db file does not exist it is created automatically.
% Re-opening an existing file is safe — createTables() is idempotent.

db = SMTDatabase('smt_pipeline.db');

% =========================================================================
%% PART A: REGISTERING A NEW EXPERIMENT SET
% =========================================================================

%% 2.  Build a TrajectoryCollection from CSV or SMD files
%
% Here we use CSV files.  The same steps work with tc.addFromSMD().

DATA_DIR  = 'path/to/csv/files';   % <-- edit
TOML_FILE = 'path/to/experiment.toml';

tc = TrajectoryCollection();
tc.ReadParams(TOML_FILE);

files = dir(fullfile(DATA_DIR, '*.csv'));
for k = 1:numel(files)
    tc.addFromFile(fullfile(files(k).folder, files(k).name), ...
                   'FileID',    num2str(k), ...
                   'Condition', 'Control');
end
tc.summary();

%% 3.  Preview what metadata will be auto-parsed from the file paths
%
% The database infers Molecule, Substrate, ImagingHz, and ExposureTime from
% the HILO directory convention and filename tokens (_e10ms_i200ms_, etc.).
% Call parsePathMetadata on any one file path to see what will be found
% before committing to the database.

if ~isempty(files)
    meta = SMTDatabase.parsePathMetadata(fullfile(files(1).folder, files(1).name));
    disp(meta)
    % → meta.molecule, meta.substrate, meta.imaging_hz,
    %   meta.exposure_time_s, meta.exp_date
end

%% 4.  Register the collection
%
% One collection groups all FOVs belonging to the same molecule / treatment
% / substrate / imaging-rate combination — across any number of acquisition
% dates.  Name it descriptively; it must be unique in the database.
%
% Treatment is the one field that cannot be inferred from file paths and
% must always be provided.  Everything else is auto-parsed (and can be
% overridden with explicit keyword arguments if the paths are non-standard).

db.registerCollection(tc, ...
    'Name',      'H2B_Control_1kPa_slow', ...
    'Treatment', 'Control');

% With all fields explicit (use when path parsing cannot be relied on):
%
%   db.registerCollection(tc, ...
%       'Name',         'H2B_Control_1kPa_slow', ...
%       'Treatment',    'Control', ...
%       'Molecule',     'H2B',  ...
%       'Substrate',    '1kPa', ...
%       'ImagingHz',    5,      ...
%       'ExposureTime', 0.010);

%% 5.  Run analyses
%
% Run whichever analyses apply.  logAnalysis() (Step 6) extracts summaries
% from the tc object, so call it right after each analysis completes.

tc.getRLDecomposition('LagTime', 4, 'String', 'H2B Control');

tc.getBayesianDiffusivity('Condition', 'Control');
% tc.computeBootstrapCI();   % optional; takes longer

%% 6.  Save the collection and register the .mat path + analysis runs
%
% tc.save() writes the full TrajectoryCollection (including results) to disk.
% db.updateCollection() stores the .mat location so loadCollection() can
% find it later.  db.logAnalysis() records a summary of each analysis run.

SAVE_PATH = 'results/H2B_Control_1kPa_slow.mat';
tc.save(SAVE_PATH);

db.updateCollection('H2B_Control_1kPa_slow', 'MatPath', SAVE_PATH);
db.logAnalysis('H2B_Control_1kPa_slow', tc, 'MSD');
db.logAnalysis('H2B_Control_1kPa_slow', tc, 'pEM');

% =========================================================================
%% PART B: BROWSING AND RELOADING ACROSS SESSIONS
% =========================================================================

%% 7.  List all registered collections
%
% With no arguments, lists every collection.
% Filter by any combination of Molecule, Treatment, Substrate, ImagingHz.

db.listCollections()
db.listCollections('Molecule', 'H2B')
db.listCollections('Molecule', 'H2B', 'Substrate', '1kPa')

%% 8.  List analysis runs
%
% Shows what has been computed and when, with key numerical summaries.

db.listAnalyses()
db.listAnalyses('AnalysisType', 'pEM')
db.listAnalyses('CollectionName', 'H2B_Control_1kPa_slow')

%% 9.  Get full details for one collection
%
% Returns a struct with:
%   .collection  – metadata row (molecule, treatment, substrate, ...)
%   .experiments – one row per FOV (date, n_nuclei, n_culled_tracks, ...)
%   .analyses    – analysis run history, newest first
%   .pEM         – decoded pEM state summaries from the most recent pEM run
%   .msd         – decoded MSD summary from the most recent MSD run

info = db.getCollection('H2B_Control_1kPa_slow');
disp(info.collection)
disp(info.experiments)
disp(info.pEM)

%% 10.  Reload a saved collection from disk
%
% Reads the .mat path stored on the collection record and loads the
% TrajectoryCollection object.  All analysis results are preserved.

tc2 = db.loadCollection('H2B_Control_1kPa_slow');
tc2.summary();

%% 11.  Raw SQL queries
%
% For anything not covered by the helper methods, use db.query() directly.

% All pEM runs with 3 diffusion states:
T = db.query('SELECT * FROM analysis_runs WHERE optimal_k = 3');
disp(T)

% Compare diffusion coefficients across substrates:
T = db.query(['SELECT c.name, c.substrate, ar.states_json ' ...
              'FROM analysis_runs ar ' ...
              'JOIN collections c ON c.id = ar.collection_id ' ...
              'WHERE c.molecule = ''H2B'' AND ar.analysis_type = ''PEM''']);
disp(T)

%% 12.  Clean up
db.close();
