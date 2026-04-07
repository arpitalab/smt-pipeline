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

tc.getMSD('LagTime', 4, 'String', 'H2B Control');

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
