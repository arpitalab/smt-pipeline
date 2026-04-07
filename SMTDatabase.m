classdef SMTDatabase < handle
    % SMTDATABASE  SQLite-backed registry for SMT pipeline experiments and analyses.
    %
    % Hierarchy
    % ---------
    %   Collection  – one per molecule / treatment / substrate / imaging rate.
    %                 Spans ALL dates and FOVs for that condition group.
    %                 Example: "H2B_E2_1kPa_slow"
    %
    %   Experiment  – one per FOV (TrajectoryWrapper).
    %                 Experiment date auto-extracted from the file path.
    %                 Many experiments → one collection (via collection_experiments).
    %
    %   Analysis run – one per completed analysis (pEM, MSD, Lifetime).
    %                  Stores key numerical summaries as JSON; full results
    %                  live in the .mat pointed to by the collection.
    %
    % Requires MATLAB R2022b+ (built-in sqlite, no Database Toolbox needed).
    %
    % -------------------------------------------------------------------------
    % Typical workflow
    % -------------------------------------------------------------------------
    %   db = SMTDatabase('smt_pipeline.db');
    %
    %   % Register once per molecule/condition group (first time or after adding data):
    %   db.registerCollection(tc, ...
    %       'Name',      'H2B_E2_1kPa_slow', ...
    %       'Treatment', 'E2');
    %   % Molecule, Substrate, ImagingHz, ExposureTime, and per-FOV dates are
    %   % auto-parsed from the file paths (HILO/{molecule}/{e-i}/{substrate}/{YYYYMMDD}/...).
    %   % Treatment is the only field that must be supplied manually.
    %
    %   % Preview what will be parsed before registering:
    %   SMTDatabase.parsePathMetadata('/path/to/MCF7_H2B_..._e10ms_i200ms_..._trajs.csv')
    %
    %   % After running analysis and saving:
    %   tc.getBayesianDiffusivity();
    %   tc.computeBootstrapCI();
    %   tc.save('results/H2B_E2_1kPa_slow.mat');
    %   db.updateCollection('H2B_E2_1kPa_slow', 'MatPath', 'results/H2B_E2_1kPa_slow.mat');
    %   db.logAnalysis('H2B_E2_1kPa_slow', tc, 'pEM');
    %   db.logAnalysis('H2B_E2_1kPa_slow', tc, 'Lifetime');
    %
    %   % Browse:
    %   db.listCollections()
    %   db.listCollections('Molecule','H2B','Substrate','1kPa')
    %   db.listAnalyses('AnalysisType','pEM')
    %   info = db.getCollection('H2B_E2_1kPa_slow');
    %       % → info.collection, info.experiments (with exp_date per FOV),
    %       %   info.analyses, info.pEM (decoded state summaries)
    %   tc   = db.loadCollection('H2B_E2_1kPa_slow');
    %   T    = db.query('SELECT * FROM analysis_runs WHERE optimal_k = 3');

    properties (Access = private)
        Conn    % sqlite connection handle
        DBPath  % path to .db file
    end

    % =========================================================================
    % PUBLIC METHODS
    % =========================================================================
    methods

        %% Constructor
        function obj = SMTDatabase(dbPath)
            % SMTDATABASE  Open (or create) an SMT experiment database.
            %
            %   db = SMTDatabase('smt_pipeline.db')
            if nargin < 1
                dbPath = 'smt_pipeline.db';
            end
            obj.DBPath = char(dbPath);
            if isfile(obj.DBPath)
                obj.Conn = sqlite(obj.DBPath, "connect");
            else
                obj.Conn = sqlite(obj.DBPath, "create");
            end
            obj.createTables();
            fprintf('SMTDatabase: connected to %s\n', obj.DBPath);
        end

        %% Destructor
        function delete(obj)
            if ~isempty(obj.Conn)
                try, close(obj.Conn); catch, end
            end
        end

        %% Explicit close
        function close(obj)
            if ~isempty(obj.Conn)
                close(obj.Conn);
                obj.Conn = [];
            end
        end

        % -----------------------------------------------------------------
        %% Register a collection
        % -----------------------------------------------------------------
        function collId = registerCollection(obj, tc, varargin)
            % REGISTERCOLLECTION  Register a TrajectoryCollection and all its FOVs.
            %
            %   Minimal call (all physical metadata auto-parsed from file paths):
            %     collId = db.registerCollection(tc, ...
            %         'Name',      'H2B_E2_1kPa_slow', ...
            %         'Treatment', 'E2')
            %
            %   Full call (override any auto-parsed field):
            %     collId = db.registerCollection(tc, ...
            %         'Name',         'H2B_E2_1kPa_slow', ...
            %         'Treatment',    'E2', ...
            %         'Molecule',     'H2B',   ...  % auto-parsed if omitted
            %         'Substrate',    '1kPa',  ...  % auto-parsed if omitted
            %         'ImagingHz',    5,       ...  % auto-parsed if omitted
            %         'ExposureTime', 0.01,    ...  % auto-parsed if omitted
            %         'MatPath',      '',      ...  % set later via updateCollection
            %         'TomlFile',     '',      ...
            %         'Notes',        '')
            %
            %   Parsing is done from the HILO directory structure and filenames:
            %     HILO/{molecule}/{e-i}/{substrate}/{YYYYMMDD}/quot_input_links/*.csv
            %   Use SMTDatabase.parsePathMetadata(path) to preview what will be found.
            %
            %   Treatment is the only field that cannot be inferred from the path
            %   and must always be provided.

            p = inputParser;
            addParameter(p, 'Name',         '',  @(x) ischar(x)||isstring(x));
            addParameter(p, 'Treatment',    '',  @(x) ischar(x)||isstring(x));
            addParameter(p, 'Molecule',     '',  @(x) ischar(x)||isstring(x));
            addParameter(p, 'Substrate',    '',  @(x) ischar(x)||isstring(x));
            addParameter(p, 'ImagingHz',    NaN, @isnumeric);
            addParameter(p, 'ExposureTime', NaN, @isnumeric);
            addParameter(p, 'MatPath',      '',  @(x) ischar(x)||isstring(x));
            addParameter(p, 'TomlFile',     '',  @(x) ischar(x)||isstring(x));
            addParameter(p, 'Notes',        '',  @(x) ischar(x)||isstring(x));
            parse(p, varargin{:});
            o = p.Results;

            name = char(o.Name);
            if isempty(name)
                error('SMTDatabase:missingName', ...
                    'registerCollection requires a ''Name'' argument.');
            end

            existing = obj.fetchOne("SELECT id FROM collections WHERE name='%s'", name);
            if ~isempty(existing)
                error('SMTDatabase:duplicateName', ...
                    ['Collection "%s" already exists.\n' ...
                     'Use updateCollection() to change metadata, or choose a new name.'], name);
            end

            nW = numel(tc.Wrappers);
            if nW == 0
                error('SMTDatabase:emptyCollection', 'TrajectoryCollection has no wrappers.');
            end

            % --- Auto-infer physical metadata from file paths ---
            % Parse metadata from every wrapper's path, then check consistency.
            % Explicit arguments always win over inferred values.
            allParsed = cell(nW, 1);
            for w = 1:nW
                fpath = char(tc.Metadata.FilePath(w));
                allParsed{w} = SMTDatabase.parsePathMetadata(fpath);
            end

            molecule     = char(o.Molecule);
            substrate    = char(o.Substrate);
            imagingHz    = o.ImagingHz;
            exposureTime = o.ExposureTime;

            if isempty(molecule)
                molecule = SMTDatabase.consensusField(allParsed, 'molecule', 'Molecule');
            end
            if isempty(substrate)
                substrate = SMTDatabase.consensusField(allParsed, 'substrate', 'Substrate');
            end
            if isnan(imagingHz)
                imagingHz = SMTDatabase.consensusNumField(allParsed, 'imaging_hz', 'ImagingHz');
            end
            if isnan(exposureTime)
                exposureTime = SMTDatabase.consensusNumField(allParsed, 'exposure_time_s', 'ExposureTime');
            end

            execute(obj.Conn, 'BEGIN TRANSACTION');
            try
                % 1. Upsert each FOV — date auto-extracted from file path
                expIds = zeros(nW, 1);
                for w = 1:nW
                    tw         = tc.Wrappers{w};
                    fid        = char(tc.Metadata.FileID(w));
                    fpath      = char(tc.Metadata.FilePath(w));
                    expIds(w)  = obj.upsertExperiment(tw, fid, fpath);
                end

                % 2. Insert collection row
                totalCulled = sum(cellfun(@(tw) tw.getNumCulledTracks(), tc.Wrappers));
                now_str     = datestr(now, 'yyyy-mm-dd HH:MM:SS');

                execute(obj.Conn, sprintf( ...
                    ['INSERT INTO collections ' ...
                     '(name,molecule,treatment,substrate,imaging_hz,exposure_time_s,' ...
                     ' n_fovs,n_culled_tracks,mat_path,toml_file,' ...
                     ' date_created,date_modified,notes) ' ...
                     'VALUES (''%s'',''%s'',''%s'',''%s'',%s,%s,' ...
                     '%d,%d,''%s'',''%s'',''%s'',''%s'',''%s'')'], ...
                    esc(name), esc(molecule), esc(o.Treatment), ...
                    esc(substrate), n2s(imagingHz), n2s(exposureTime), ...
                    nW, totalCulled, esc(o.MatPath), esc(o.TomlFile), ...
                    now_str, now_str, esc(o.Notes)));

                res2   = fetch(obj.Conn, 'SELECT last_insert_rowid() AS id');
                collId = res2.id(1);

                % 3. Link FOVs to collection, preserving condition tags
                for w = 1:nW
                    condTag = char(tc.Metadata.Condition(w));
                    execute(obj.Conn, sprintf( ...
                        ['INSERT OR IGNORE INTO collection_experiments ' ...
                         '(collection_id,experiment_id,condition_tag) VALUES (%d,%d,''%s'')'], ...
                        collId, expIds(w), esc(condTag)));
                end

                execute(obj.Conn, 'COMMIT');

            catch ME
                execute(obj.Conn, 'ROLLBACK');
                rethrow(ME);
            end

            % Report date range across FOVs
            dateRange = obj.collectionDateRange(collId);
            fprintf('Registered "%s" (id=%d): %d FOVs, %d culled tracks, dates: %s\n', ...
                name, collId, nW, totalCulled, dateRange);
        end

        % -----------------------------------------------------------------
        %% Update collection metadata
        % -----------------------------------------------------------------
        function updateCollection(obj, name, varargin)
            % UPDATECOLLECTION  Update fields on an existing collection.
            %
            %   db.updateCollection('H2B_E2_1kPa_slow', ...
            %       'MatPath',      'results/H2B_E2_1kPa_slow.mat', ...
            %       'ExposureTime', 0.01, ...
            %       'Notes',        're-analyzed March 2024');
            %
            %   Updatable: MatPath, Notes, Molecule, Treatment, Substrate,
            %              ImagingHz, ExposureTime, TomlFile

            p = inputParser;
            addParameter(p, 'MatPath',      [], @(x) ischar(x)||isstring(x));
            addParameter(p, 'Notes',        [], @(x) ischar(x)||isstring(x));
            addParameter(p, 'Molecule',     [], @(x) ischar(x)||isstring(x));
            addParameter(p, 'Treatment',    [], @(x) ischar(x)||isstring(x));
            addParameter(p, 'Substrate',    [], @(x) ischar(x)||isstring(x));
            addParameter(p, 'TomlFile',     [], @(x) ischar(x)||isstring(x));
            addParameter(p, 'ImagingHz',    [], @isnumeric);
            addParameter(p, 'ExposureTime', [], @isnumeric);
            parse(p, varargin{:});
            o = p.Results;

            strFields = {'MatPath','mat_path'; 'Notes','notes'; 'Molecule','molecule'; ...
                         'Treatment','treatment'; 'Substrate','substrate'; 'TomlFile','toml_file'};
            sets = {};
            for f = 1:size(strFields,1)
                v = o.(strFields{f,1});
                if ~isempty(v)
                    sets{end+1} = sprintf("%s='%s'", strFields{f,2}, esc(v)); %#ok<AGROW>
                end
            end
            for nf = {'ImagingHz','imaging_hz'; 'ExposureTime','exposure_time_s'}'
                v = o.(nf{1});
                if ~isempty(v)
                    sets{end+1} = sprintf('%s=%s', nf{2}, n2s(v)); %#ok<AGROW>
                end
            end
            if isempty(sets), return; end

            sets{end+1} = sprintf("date_modified='%s'", datestr(now,'yyyy-mm-dd HH:MM:SS'));
            execute(obj.Conn, sprintf("UPDATE collections SET %s WHERE name='%s'", ...
                strjoin(sets, ', '), esc(name)));
            fprintf('Updated collection "%s".\n', char(name));
        end

        % -----------------------------------------------------------------
        %% Log an analysis run (manual — call after verifying results)
        % -----------------------------------------------------------------
        function runId = logAnalysis(obj, collectionName, tc, analysisType, varargin)
            % LOGANALYSIS  Record a completed analysis in the database.
            %
            %   db.logAnalysis('H2B_E2_1kPa_slow', tc, 'pEM')
            %   db.logAnalysis('H2B_E2_1kPa_slow', tc, 'pEM', ...
            %       'ConditionFilter', 'E2_treated', ...
            %       'Notes',           'K=3, bootstrap 500 resamples')
            %   db.logAnalysis('H2B_E2_1kPa_slow', tc, 'MSD')
            %   db.logAnalysis('H2B_E2_1kPa_slow', tc, 'Lifetime')
            %
            %   The mat_path stored on the collection is automatically included.
            %   analysisType : 'pEM' | 'MSD' | 'Lifetime'

            p = inputParser;
            addParameter(p, 'ConditionFilter', '', @(x) ischar(x)||isstring(x));
            addParameter(p, 'Notes',           '', @(x) ischar(x)||isstring(x));
            parse(p, varargin{:});
            o = p.Results;

            res = obj.fetchOne("SELECT id, mat_path FROM collections WHERE name='%s'", ...
                               collectionName);
            if isempty(res)
                error('SMTDatabase:notFound', ...
                    'No collection named "%s". Call registerCollection first.', ...
                    char(collectionName));
            end
            collId  = res.id;
            matPath = char(res.mat_path);

            atype = upper(strtrim(char(analysisType)));
            row   = struct();
            row.collection_id    = collId;
            row.analysis_type    = atype;
            row.condition_filter = char(o.ConditionFilter);
            row.date_run         = datestr(now, 'yyyy-mm-dd HH:MM:SS');
            row.mat_path         = matPath;
            row.notes            = char(o.Notes);

            switch atype
                case 'PEM',      row = SMTDatabase.extractPEMSummary(row, tc);
                case 'MSD',      row = SMTDatabase.extractMSDSummary(row, tc);
                case 'LIFETIME', row = SMTDatabase.extractLifetimeSummary(row, tc);
                otherwise
                    warning('SMTDatabase:unknownType', ...
                        'Unknown analysisType "%s". Logging metadata only.', atype);
            end

            flds     = fieldnames(row);
            vals     = struct2cell(row);
            colStr   = strjoin(flds, ', ');
            valParts = cell(numel(vals), 1);
            for i = 1:numel(vals)
                v = vals{i};
                if isnumeric(v) && isscalar(v)
                    valParts{i} = n2s(v);
                else
                    valParts{i} = sprintf("'%s'", esc(char(v)));
                end
            end
            execute(obj.Conn, sprintf('INSERT INTO analysis_runs (%s) VALUES (%s)', ...
                colStr, strjoin(valParts, ', ')));

            res2  = fetch(obj.Conn, 'SELECT last_insert_rowid() AS id');
            runId = res2.id(1);
            fprintf('Logged %s for "%s" (run id=%d).\n', atype, char(collectionName), runId);
        end

        % -----------------------------------------------------------------
        %% List collections
        % -----------------------------------------------------------------
        function T = listCollections(obj, varargin)
            % LISTCOLLECTIONS  Show registered collections, optionally filtered.
            %
            %   db.listCollections()
            %   T = db.listCollections('Molecule','H2B','Substrate','1kPa')
            %   T = db.listCollections('ImagingHz', 100)

            p = inputParser;
            addParameter(p, 'Molecule',  '', @(x) ischar(x)||isstring(x));
            addParameter(p, 'Treatment', '', @(x) ischar(x)||isstring(x));
            addParameter(p, 'Substrate', '', @(x) ischar(x)||isstring(x));
            addParameter(p, 'ImagingHz', NaN, @isnumeric);
            parse(p, varargin{:});
            o = p.Results;

            parts = {};
            if ~isempty(char(o.Molecule)),  parts{end+1} = sprintf("c.molecule='%s'",  esc(o.Molecule));  end
            if ~isempty(char(o.Treatment)), parts{end+1} = sprintf("c.treatment='%s'", esc(o.Treatment)); end
            if ~isempty(char(o.Substrate)), parts{end+1} = sprintf("c.substrate='%s'", esc(o.Substrate)); end
            if ~isnan(o.ImagingHz),         parts{end+1} = sprintf('c.imaging_hz=%s',  n2s(o.ImagingHz)); end
            w = SMTDatabase.partsToWhere(parts);

            % Date range comes from linked experiments, not from the collection row
            sql = ['SELECT c.id, c.name, c.molecule, c.treatment, c.substrate, ' ...
                   'c.imaging_hz, c.exposure_time_s, c.n_fovs, c.n_culled_tracks, ' ...
                   'MIN(e.exp_date) AS first_date, MAX(e.exp_date) AS last_date, ' ...
                   'COUNT(DISTINCT e.exp_date) AS n_dates, ' ...
                   'COUNT(DISTINCT ar.id) AS n_analyses, ' ...
                   'c.mat_path, c.date_created, c.notes ' ...
                   'FROM collections c ' ...
                   'LEFT JOIN collection_experiments ce ON ce.collection_id = c.id ' ...
                   'LEFT JOIN experiments e  ON e.id  = ce.experiment_id ' ...
                   'LEFT JOIN analysis_runs ar ON ar.collection_id = c.id ' ...
                   w ' GROUP BY c.id ORDER BY c.molecule, c.treatment, c.substrate, c.imaging_hz'];
            T = fetch(obj.Conn, sql);
            if nargout == 0
                if isempty(T) || height(T) == 0
                    fprintf('No collections registered yet.\n');
                else
                    disp(T(:, {'id','name','molecule','treatment','substrate', ...
                               'imaging_hz','n_fovs','n_culled_tracks', ...
                               'first_date','last_date','n_dates','n_analyses'}));
                end
            end
        end

        % -----------------------------------------------------------------
        %% List analysis runs
        % -----------------------------------------------------------------
        function T = listAnalyses(obj, varargin)
            % LISTANALYSES  Show analysis runs, optionally filtered.
            %
            %   db.listAnalyses()
            %   T = db.listAnalyses('Molecule','H2B','AnalysisType','pEM')
            %   T = db.listAnalyses('CollectionName','H2B_E2_1kPa_slow')

            p = inputParser;
            addParameter(p, 'Molecule',       '', @(x) ischar(x)||isstring(x));
            addParameter(p, 'Treatment',      '', @(x) ischar(x)||isstring(x));
            addParameter(p, 'Substrate',      '', @(x) ischar(x)||isstring(x));
            addParameter(p, 'AnalysisType',   '', @(x) ischar(x)||isstring(x));
            addParameter(p, 'CollectionName', '', @(x) ischar(x)||isstring(x));
            parse(p, varargin{:});
            o = p.Results;

            parts = {};
            if ~isempty(char(o.Molecule)),       parts{end+1} = sprintf("c.molecule='%s'",       esc(o.Molecule));    end
            if ~isempty(char(o.Treatment)),      parts{end+1} = sprintf("c.treatment='%s'",      esc(o.Treatment));   end
            if ~isempty(char(o.Substrate)),      parts{end+1} = sprintf("c.substrate='%s'",      esc(o.Substrate));   end
            if ~isempty(char(o.AnalysisType)),   parts{end+1} = sprintf("ar.analysis_type='%s'", esc(upper(char(o.AnalysisType)))); end
            if ~isempty(char(o.CollectionName)), parts{end+1} = sprintf("c.name='%s'",           esc(o.CollectionName)); end
            w = SMTDatabase.partsToWhere(parts);

            sql = ['SELECT ar.id, c.name AS collection, c.molecule, c.treatment, ' ...
                   'c.substrate, ar.analysis_type, ar.condition_filter, ar.date_run, ' ...
                   'ar.optimal_k, ar.n_raw_tracks, ar.bootstrap_done, ar.mat_path, ar.notes ' ...
                   'FROM analysis_runs ar ' ...
                   'JOIN collections c ON c.id = ar.collection_id ' ...
                   w ' ORDER BY ar.date_run DESC'];
            T = fetch(obj.Conn, sql);
            if nargout == 0
                if isempty(T) || height(T) == 0
                    fprintf('No analysis runs found.\n');
                else
                    disp(T);
                end
            end
        end

        % -----------------------------------------------------------------
        %% Get all info for a named collection
        % -----------------------------------------------------------------
        function info = getCollection(obj, name)
            % GETCOLLECTION  Return all info for a named collection as a struct.
            %
            %   info = db.getCollection('H2B_E2_1kPa_slow')
            %
            %   info.collection  – 1-row table: name, molecule, treatment, ...
            %   info.experiments – table: one row per FOV, includes exp_date,
            %                      n_nuclei, n_culled_tracks, condition_tag
            %   info.analyses    – table: analysis runs, newest first
            %   info.pEM         – decoded pEM state summary (last run)
            %   info.msd         – decoded MSD summary (last run)
            %   info.lifetime    – decoded Lifetime summary (last run)

            T = fetch(obj.Conn, ...
                sprintf("SELECT * FROM collections WHERE name='%s'", esc(name)));
            if isempty(T) || height(T) == 0
                error('SMTDatabase:notFound', 'No collection named "%s".', char(name));
            end
            info.collection = T;
            collId = T.id(1);

            info.experiments = fetch(obj.Conn, sprintf( ...
                ['SELECT e.file_id, e.exp_date, e.source_type, ' ...
                 'e.pixel_size_um, e.frame_interval_s, ' ...
                 'e.n_nuclei, e.n_raw_tracks, e.n_culled_tracks, ' ...
                 'ce.condition_tag, e.file_path ' ...
                 'FROM experiments e ' ...
                 'JOIN collection_experiments ce ON ce.experiment_id = e.id ' ...
                 'WHERE ce.collection_id = %d ORDER BY e.exp_date, e.file_id'], collId));

            info.analyses = fetch(obj.Conn, sprintf( ...
                ['SELECT * FROM analysis_runs WHERE collection_id=%d ' ...
                 'ORDER BY date_run DESC'], collId));

            % Decode JSON from the most recent run of each type
            info.pEM      = SMTDatabase.decodeLatestJSON(info.analyses, 'PEM',      'states_json');
            info.msd      = SMTDatabase.decodeLatestJSON(info.analyses, 'MSD',      'msd_json');
            info.lifetime = SMTDatabase.decodeLatestJSON(info.analyses, 'LIFETIME', 'lifetime_json');
        end

        % -----------------------------------------------------------------
        %% Load the saved .mat for a named collection
        % -----------------------------------------------------------------
        function tc = loadCollection(obj, name)
            % LOADCOLLECTION  Load the saved .mat for a named collection.
            %
            %   tc = db.loadCollection('H2B_E2_1kPa_slow')

            res = obj.fetchOne("SELECT mat_path FROM collections WHERE name='%s'", name);
            if isempty(res)
                error('SMTDatabase:notFound', 'No collection named "%s".', char(name));
            end
            matPath = char(res.mat_path);
            if isempty(matPath)
                error('SMTDatabase:noMatPath', ...
                    ['No mat_path set for "%s".\n' ...
                     'Run tc.save() then db.updateCollection(''%s'',''MatPath'',path).'], ...
                    char(name), char(name));
            end
            if ~isfile(matPath)
                error('SMTDatabase:fileNotFound', '.mat not found: %s', matPath);
            end
            data = load(matPath, 'tc');
            tc   = data.tc;
            fprintf('Loaded "%s" from %s\n', char(name), matPath);
        end

        % -----------------------------------------------------------------
        %% Raw SQL passthrough
        % -----------------------------------------------------------------
        function T = query(obj, sql)
            % QUERY  Execute arbitrary SQL and return results as a table.
            %
            %   T = db.query('SELECT * FROM collections WHERE substrate = ''Glass''')
            %   T = db.query(['SELECT c.name, ar.optimal_k, ar.states_json ' ...
            %                 'FROM analysis_runs ar JOIN collections c ON c.id=ar.collection_id ' ...
            %                 'WHERE c.molecule=''H2B'' AND ar.analysis_type=''PEM'''])
            T = fetch(obj.Conn, sql);
        end

        % -----------------------------------------------------------------
        %% Delete a collection record (does not touch the .mat on disk)
        % -----------------------------------------------------------------
        function deleteCollection(obj, name)
            % DELETECOLLECTION  Remove a collection and its analysis runs from the DB.
            %
            %   db.deleteCollection('H2B_E2_1kPa_slow')
            %
            %   Removes: the collection row, its analysis_runs, and the FOV links.
            %   Does NOT remove experiment rows (shared across collections) or the .mat.

            res = obj.fetchOne("SELECT id FROM collections WHERE name='%s'", name);
            if isempty(res)
                warning('SMTDatabase:notFound', 'No collection named "%s".', char(name));
                return;
            end
            collId = res.id;

            execute(obj.Conn, 'BEGIN TRANSACTION');
            try
                execute(obj.Conn, sprintf('DELETE FROM analysis_runs          WHERE collection_id=%d', collId));
                execute(obj.Conn, sprintf('DELETE FROM collection_experiments WHERE collection_id=%d', collId));
                execute(obj.Conn, sprintf('DELETE FROM collections            WHERE id=%d',            collId));
                execute(obj.Conn, 'COMMIT');
            catch ME
                execute(obj.Conn, 'ROLLBACK');
                rethrow(ME);
            end
            fprintf('Deleted collection "%s" and its analysis runs from the database.\n', char(name));
        end

    end % public methods

    % =========================================================================
    % PRIVATE METHODS
    % =========================================================================
    methods (Access = private)

        %% Create tables (idempotent — safe to call on every open)
        function createTables(obj)
            execute(obj.Conn, [ ...
                'CREATE TABLE IF NOT EXISTS experiments (' ...
                '  id               INTEGER PRIMARY KEY AUTOINCREMENT,' ...
                '  file_id          TEXT,' ...
                '  file_path        TEXT UNIQUE,' ...
                '  source_type      TEXT,' ...
                '  pixel_size_um    REAL,' ...
                '  frame_interval_s REAL,' ...
                '  n_nuclei         INTEGER,' ...
                '  n_raw_tracks     INTEGER,' ...
                '  n_culled_tracks  INTEGER,' ...
                '  exp_date         TEXT,' ...    % auto-extracted from file path
                '  date_registered  TEXT,' ...
                '  notes            TEXT' ...
                ')']);

            execute(obj.Conn, [ ...
                'CREATE TABLE IF NOT EXISTS collections (' ...
                '  id               INTEGER PRIMARY KEY AUTOINCREMENT,' ...
                '  name             TEXT UNIQUE NOT NULL,' ...
                '  molecule         TEXT,' ...
                '  treatment        TEXT,' ...
                '  substrate        TEXT,' ...
                '  imaging_hz       REAL,' ...
                '  exposure_time_s  REAL,' ...
                '  n_fovs           INTEGER,' ...
                '  n_culled_tracks  INTEGER,' ...
                '  mat_path         TEXT,' ...
                '  toml_file        TEXT,' ...
                '  date_created     TEXT,' ...
                '  date_modified    TEXT,' ...
                '  notes            TEXT' ...
                ')']);
            % Note: no exp_date on collections — dates live on experiments
            % and are aggregated (first_date / last_date) in listCollections.

            execute(obj.Conn, [ ...
                'CREATE TABLE IF NOT EXISTS collection_experiments (' ...
                '  collection_id INTEGER REFERENCES collections(id),' ...
                '  experiment_id INTEGER REFERENCES experiments(id),' ...
                '  condition_tag TEXT,' ...   % the tc.Metadata.Condition value for this FOV
                '  PRIMARY KEY (collection_id, experiment_id)' ...
                ')']);

            execute(obj.Conn, [ ...
                'CREATE TABLE IF NOT EXISTS analysis_runs (' ...
                '  id               INTEGER PRIMARY KEY AUTOINCREMENT,' ...
                '  collection_id    INTEGER REFERENCES collections(id),' ...
                '  analysis_type    TEXT,' ...           % PEM | MSD | LIFETIME
                '  condition_filter TEXT,' ...           % empty = full collection
                '  date_run         TEXT,' ...
                '  mat_path         TEXT,' ...
                '  parameters_json  TEXT,' ...           % pEM: dt, R, splitLength, lambda, maxStates
                '  optimal_k        INTEGER,' ...        % pEM: number of states
                '  log_likelihood   REAL,' ...           % pEM: optimal log-likelihood
                '  n_raw_tracks     INTEGER,' ...        % pEM: parent tracks fed in
                '  n_split_tracks   INTEGER,' ...        % pEM: split-track count
                '  states_json      TEXT,' ...           % pEM: [{state,D,sigma,proportion,D_ci_lo,...}]
                '  bootstrap_done   INTEGER DEFAULT 0,' ...
                '  bootstrap_nboot  INTEGER,' ...
                '  msd_json         TEXT,' ...           % MSD: {M:[...]}
                '  lifetime_json    TEXT,' ...           % Lifetime: {x,f,flo,fup}
                '  notes            TEXT' ...
                ')']);
        end

        %% Upsert a single FOV into experiments; return its row id
        function expId = upsertExperiment(obj, tw, fileId, filePath)
            if isempty(filePath)
                filePath = char(tw.FileName);
            end

            res     = obj.fetchOne("SELECT id FROM experiments WHERE file_path='%s'", filePath);
            nNuclei = numel(tw.ROIs);
            nRaw    = tw.getNumRawTracks();
            nCulled = tw.getNumCulledTracks();
            ps      = tw.getPixelSize();
            fi      = tw.getFrameInterval();
            st      = tw.getSourceType();
            expDate = SMTDatabase.extractDateFromPath(filePath);

            if ~isempty(res)
                % Refresh counts in case tracks were re-processed
                expId = res.id;
                execute(obj.Conn, sprintf( ...
                    ['UPDATE experiments SET ' ...
                     'n_nuclei=%d, n_raw_tracks=%d, n_culled_tracks=%d, ' ...
                     'pixel_size_um=%s, frame_interval_s=%s ' ...
                     'WHERE id=%d'], ...
                    nNuclei, nRaw, nCulled, n2s(ps), n2s(fi), expId));
            else
                execute(obj.Conn, sprintf( ...
                    ['INSERT INTO experiments ' ...
                     '(file_id,file_path,source_type,pixel_size_um,frame_interval_s,' ...
                     ' n_nuclei,n_raw_tracks,n_culled_tracks,exp_date,date_registered) ' ...
                     'VALUES (''%s'',''%s'',''%s'',%s,%s,%d,%d,%d,''%s'',''%s'')'], ...
                    esc(fileId), esc(filePath), esc(st), n2s(ps), n2s(fi), ...
                    nNuclei, nRaw, nCulled, ...
                    expDate, datestr(now,'yyyy-mm-dd HH:MM:SS')));
                res2  = fetch(obj.Conn, 'SELECT last_insert_rowid() AS id');
                expId = res2.id(1);
            end
        end

        %% Build date-range string for display
        function s = collectionDateRange(obj, collId)
            res = fetch(obj.Conn, sprintf( ...
                ['SELECT MIN(e.exp_date) AS d0, MAX(e.exp_date) AS d1, ' ...
                 'COUNT(DISTINCT e.exp_date) AS nd ' ...
                 'FROM collection_experiments ce ' ...
                 'JOIN experiments e ON e.id = ce.experiment_id ' ...
                 'WHERE ce.collection_id = %d'], collId));
            if isempty(res) || height(res) == 0 || isempty(char(res.d0(1)))
                s = 'unknown (no date in filenames)';
            elseif strcmp(char(res.d0(1)), char(res.d1(1)))
                s = sprintf('%s (%d replicates)', char(res.d0(1)), res.nd(1));
            else
                s = sprintf('%s → %s (%d dates)', char(res.d0(1)), char(res.d1(1)), res.nd(1));
            end
        end

        %% Fetch one row as a scalar struct (empty if not found)
        function s = fetchOne(obj, fmt, varargin)
            args = cellfun(@esc, varargin, 'UniformOutput', false);
            sql  = sprintf(fmt, args{:});
            T    = fetch(obj.Conn, sql);
            if isempty(T) || height(T) == 0
                s = [];
            else
                s = table2struct(T(1,:));
            end
        end

    end % private methods

    % =========================================================================
    % STATIC PRIVATE HELPERS
    % =========================================================================
    methods (Static)  % public static — callable as SMTDatabase.parsePathMetadata(...)

        function meta = parsePathMetadata(filepath)
            % PARSEPATHMETADATA  Extract experimental metadata from a HILO file path.
            %
            %   meta = SMTDatabase.parsePathMetadata(filepath)
            %
            %   Parses the HILO directory convention:
            %     .../HILO/{molecule}/{eXms-iYms}/{substrate}/{YYYYMMDD}/quot_input_links/
            %             MCF7_{molecule}_{substrate}_..._e{E}ms_i{I}ms_..._trajs.csv
            %
            %   Returns a struct with fields:
            %     .molecule       – 'H2B' | 'ER' | '' (if not found)
            %     .substrate      – '1kPa' | '12kPa' | '100kPa' | 'Glass' | ''
            %     .imaging_hz     – frame rate in Hz  (e.g. 5 or 100), or NaN
            %     .exposure_time_s – camera exposure in seconds (e.g. 0.010), or NaN
            %     .exp_date       – 'YYYY-MM-DD' | ''
            %
            %   Example:
            %     meta = SMTDatabase.parsePathMetadata( ...
            %         '/HILO/H2B/10-200/1kPa/20241015/quot_input_links/MCF7_H2B_1kPa_..._e10ms_i200ms_..._trajs.csv')
            %     % → molecule='H2B', substrate='1kPa', imaging_hz=5, exposure_time_s=0.010, exp_date='2024-10-15'

            fp = char(filepath);

            meta.molecule        = '';
            meta.substrate       = '';
            meta.imaging_hz      = NaN;
            meta.exposure_time_s = NaN;
            meta.exp_date        = '';

            % Molecule: directory segment or filename token
            tok = regexp(fp, '[/\\](H2B|ER)[/\\]', 'tokens', 'once');
            if isempty(tok)
                tok = regexp(fp, '_(H2B|ER)_', 'tokens', 'once');
            end
            if ~isempty(tok), meta.molecule = tok{1}; end

            % Substrate: directory segment or filename token
            tok = regexp(fp, '[/\\](100kPa|12kPa|1kPa|Glass)[/\\]', 'tokens', 'once');
            if isempty(tok)
                tok = regexp(fp, '_(100kPa|12kPa|1kPa|Glass)_', 'tokens', 'once');
            end
            if ~isempty(tok), meta.substrate = tok{1}; end

            % Exposure and interval: from filename _eXms_iYms_ or directory XX-YY
            tok = regexp(fp, '_e(\d+)ms_i(\d+)ms_', 'tokens', 'once');
            if isempty(tok)
                % Fall back to directory name like "10-200" or "10-10"
                tok = regexp(fp, '[/\\](\d+)-(\d+)[/\\]', 'tokens', 'once');
            end
            if ~isempty(tok)
                eMs = str2double(tok{1});
                iMs = str2double(tok{2});
                meta.exposure_time_s = eMs / 1000;
                meta.imaging_hz      = 1000 / iMs;
            end

            % Date: YYYYMMDD directory segment (preferred) or elsewhere in path
            tok = regexp(fp, '[/\\](20\d{6})[/\\]', 'tokens', 'once');
            if ~isempty(tok)
                s = tok{1};
                meta.exp_date = sprintf('%s-%s-%s', s(1:4), s(5:6), s(7:8));
            else
                meta.exp_date = SMTDatabase.extractDateFromPath(fp);
            end
        end

    end % public static methods

    methods (Static, Access = private)

        %% Extract YYYY-MM-DD from a file path (generic fallback)
        function d = extractDateFromPath(filepath)
            d = '';
            % Compact 8-digit date (e.g. 20240315 in a path segment or filename)
            tok = regexp(filepath, '(?<!\d)(20\d{2})(0[1-9]|1[0-2])(0[1-9]|[12]\d|3[01])(?!\d)', 'tokens', 'once');
            if ~isempty(tok)
                d = sprintf('%s-%s-%s', tok{1}, tok{2}, tok{3});
                return;
            end
            % Separated date (e.g. 2024-03-15 or 2024_03_15)
            tok = regexp(filepath, '(20\d{2})[-_](0[1-9]|1[0-2])[-_](0[1-9]|[12]\d|3[01])', 'tokens', 'once');
            if ~isempty(tok)
                d = sprintf('%s-%s-%s', tok{1}, tok{2}, tok{3});
            end
        end

        %% Return consensus string value from parsed metadata; warn if inconsistent
        function val = consensusField(parsedList, field, label)
            vals = cellfun(@(m) m.(field), parsedList, 'UniformOutput', false);
            vals = vals(~cellfun(@isempty, vals));
            uvals = unique(vals);
            if isempty(uvals)
                val = '';
                warning('SMTDatabase:noMetadata', ...
                    'Could not auto-parse %s from file paths. Pass it explicitly.', label);
            elseif numel(uvals) > 1
                val = uvals{1};
                warning('SMTDatabase:inconsistentMetadata', ...
                    '%s is inconsistent across FOVs: {%s}. Using "%s". Pass it explicitly to override.', ...
                    label, strjoin(uvals,', '), val);
            else
                val = uvals{1};
            end
        end

        %% Return consensus numeric value from parsed metadata; warn if inconsistent
        function val = consensusNumField(parsedList, field, label)
            vals = cellfun(@(m) m.(field), parsedList);
            vals = vals(~isnan(vals));
            uvals = unique(vals);
            if isempty(uvals)
                val = NaN;
                warning('SMTDatabase:noMetadata', ...
                    'Could not auto-parse %s from file paths. Pass it explicitly.', label);
            elseif numel(uvals) > 1
                val = uvals(1);
                warning('SMTDatabase:inconsistentMetadata', ...
                    '%s is inconsistent across FOVs. Using %.4g. Pass it explicitly to override.', ...
                    label, val);
            else
                val = uvals(1);
            end
        end

        %% Extract pEM summary into row struct
        function row = extractPEMSummary(row, tc)
            if ~isfield(tc.pEMResults, 'pEMTable') || isempty(tc.pEMResults.pEMTable)
                warning('SMTDatabase:noPEM', 'pEM not computed on this collection.');
                return;
            end
            tbl = tc.pEMResults.pEMTable;
            K   = tbl.optimalSize(1);
            row.optimal_k      = K;
            row.log_likelihood = tbl.optimalL(1);
            row.n_raw_tracks   = tbl.numRawTracks{1};
            row.n_split_tracks = tbl.numSplitTracks(1);

            % Per-state summary (sorted by D, same order as proportions table)
            prop = tc.pEMResults.proportions;
            hasCi = isfield(tc.pEMResults,'ciresults') && ~isempty(tc.pEMResults.ciresults);
            states = struct( ...
                'state',      num2cell(prop.State), ...
                'D',          num2cell(prop.OptimalD), ...
                'sigma',      num2cell(prop.OptimalS), ...
                'proportion', num2cell(prop.Proportion));
            if hasCi
                ci = tc.pEMResults.ciresults;
                for k = 1:K
                    states(k).D_ci_lo = ci.D.ci_lo(k);
                    states(k).D_ci_hi = ci.D.ci_hi(k);
                    states(k).P_ci_lo = ci.P.ci_lo(k);
                    states(k).P_ci_hi = ci.P.ci_hi(k);
                end
                row.bootstrap_done  = 1;
                row.bootstrap_nboot = ci.nBoot;
            end
            row.states_json = jsonencode(states);

            if isfield(tc.pEMParams, 'trackInfo')
                ti   = tc.pEMParams.trackInfo;
                pars = struct('dt',ti.dt,'R',ti.R,'splitLength',ti.splitLength,'lambda',ti.lambda);
                if isfield(tc.pEMParams,'params')
                    pars.minStates = tc.pEMParams.params.minStates;
                    pars.maxStates = tc.pEMParams.params.maxStates;
                end
                row.parameters_json = jsonencode(pars);
            end
        end

        %% Extract MSD summary into row struct
        function row = extractMSDSummary(row, tc)
            if ~isfield(tc.MSDResults,'M') || isempty(tc.MSDResults.M)
                warning('SMTDatabase:noMSD', 'MSD not computed on this collection.');
                return;
            end
            row.msd_json = jsonencode(struct('M', tc.MSDResults.M(:)'));
        end

        %% Extract Lifetime summary into row struct
        function row = extractLifetimeSummary(row, tc)
            if ~isfield(tc.LifetimeResults,'x') || isempty(tc.LifetimeResults.x)
                warning('SMTDatabase:noLifetime', 'Lifetime not computed on this collection.');
                return;
            end
            lr = tc.LifetimeResults;
            row.lifetime_json = jsonencode(struct( ...
                'x', lr.x(:)', 'f', lr.f(:)', 'flo', lr.flo(:)', 'fup', lr.fup(:)'));
        end

        %% Decode JSON from the most recent run of a given analysis type
        function out = decodeLatestJSON(analyses, atype, jsonField)
            out = [];
            if isempty(analyses) || height(analyses) == 0, return; end
            idx = find(strcmp(string(analyses.analysis_type), atype), 1, 'first');
            if isempty(idx), return; end
            jsonStr = char(analyses.(jsonField)(idx));
            if ~isempty(jsonStr)
                try, out = jsondecode(jsonStr); catch, out = jsonStr; end
            end
        end

        %% Build WHERE clause from a cell array of condition strings
        function w = partsToWhere(parts)
            if isempty(parts), w = ''; else, w = ['WHERE ' strjoin(parts,' AND ')]; end
        end

    end % static private methods

end % classdef

% =============================================================================
% File-level helpers (not class members — used by module-level SQL builders)
% =============================================================================

function s = esc(v)
    % Escape single quotes for SQLite string literals.
    s = strrep(char(v), '''', '''''');
end

function s = n2s(v)
    % Convert numeric scalar to SQL literal; NaN/empty → NULL.
    if isempty(v) || (isnumeric(v) && isscalar(v) && isnan(v))
        s = 'NULL';
    else
        s = num2str(double(v), '%.10g');
    end
end
