classdef TrajectoryCollection < handle
    % TRAJECTORYCOLLECTION - Aggregates multiple TrajectoryWrapper instances for
    %                        multi-file / multi-condition SMT analysis.
    %
    % Supports two input paths:
    %   CSV path  : tc.addFromFile('locs.csv', ...)
    %   SMD path  : tc.addFromSMD(smd, ...)
    %
    % Both paths produce identical TrajectoryWrapper objects that flow into
    % the same MSD / Bayesian / lifetime analysis methods.
    %
    % Example (mixed sources):
    %   tc = TrajectoryCollection();
    %   tc.ReadParams('experiment.toml');
    %   tc.addFromFile('s001.csv', 'FileID','1', 'Condition','Control');
    %   tc.addFromSMD(smd,         'FileID','2', 'Condition','Treated');
    %   tc.getBayesianDiffusivity('Condition','Control');
    %   tc.summary();

    properties (SetAccess = private)
        Wrappers           = {}           % Cell array of TrajectoryWrapper instances
        Metadata                          % Table: FileID, FilePath, Condition
        AllRawTracks       = {}           % Lazy: aggregated raw tracks
        AllCulledTracks    = {}           % Lazy: aggregated culled tracks
        AllRelativeTracks  = {}           % Lazy: aggregated relative tracks
        RLResults          = struct()
        MSDResults         = struct()
        pEMResults         = struct()
        Params             = struct()
        pEMParams          = struct()
        LifetimeResults    = struct()
        FBMResults         = struct()    % Population fBM MLE (K, alpha, sigma)
        FBMAlphaResults    = struct()    % Per-track (K, alpha) distribution
    end

    properties (Access = private)
        IsAllRawCollected      = false
        IsAllCulledCollected   = false
        IsAllRelativeCollected = false
        IsRLComputed           = false
        IsMSDComputed          = false
        IspEMComputed          = false
        IsBootstrapComputed    = false
        IsLifetimeComputed     = false
        IsFBMComputed          = false
        IsFBMAlphaComputed     = false
        pEMBootstrapInputs     = struct()  % deltaX, rawResults, trackInfo, params
    end

    methods
        %% Constructor
        function obj = TrajectoryCollection(varargin)
            p = inputParser;
            addParameter(p, 'Condition',  '');
            addParameter(p, 'Parameters', struct());
            parse(p, varargin{:});

            obj.Metadata = table('Size', [0 3], ...
                'VariableTypes', {'string','string','string'}, ...
                'VariableNames', {'FileID','FilePath','Condition'});

            if ~isempty(p.Results.Condition)
                obj.Metadata.Condition(:) = p.Results.Condition;
            end
            if ~isempty(fieldnames(p.Results.Parameters))
                obj.addGlobalParameters(p.Results.Parameters);
            end
        end

        %% Add from CSV file
        function obj = addFromFile(obj, csvFile, varargin)
            % ADDFROMFILE  Create a TrajectoryWrapper from a CSV file, cull, and add.
            %
            %   tc.addFromFile(csvFile)
            %   tc.addFromFile(csvFile, 'Condition', 'Control', ...
            %                           'PixelSize', 0.165, ...
            %                           'FrameInterval', 0.02, ...
            %                           'Parameters', params)
            %   tc.addFromFile(csvFile, 'ParamsFile', 'experiment.toml', ...)
            %
            %   PixelSize and FrameInterval can be supplied directly or read
            %   from a TOML file via 'ParamsFile'.  Direct values take
            %   precedence over TOML values.  PixelSize must be known before
            %   readData() runs (it scales x,y from pixels to µm on load).

            p = inputParser;
            addParameter(p, 'FileID',        '');
            addParameter(p, 'Condition',      '');
            addParameter(p, 'Parameters',     struct());
            addParameter(p, 'ParamsFile',     '');       % path to TOML file
            addParameter(p, 'PixelSize',      NaN);      % µm/px; NaN = use TOML or TW default
            addParameter(p, 'FrameInterval',  NaN);      % s;     NaN = use TOML or leave unset
            parse(p, varargin{:});

            fileID = p.Results.FileID;
            if isempty(fileID)
                [~, fileID, ~] = fileparts(csvFile);
            end

            tw = TrajectoryWrapper('Parameters', p.Results.Parameters);

            % Load TOML first (lowest precedence)
            if ~isempty(p.Results.ParamsFile)
                tw.readParams(p.Results.ParamsFile);
            end

            % Direct arguments override TOML values
            if ~isnan(p.Results.PixelSize)
                tw.setPixelSize(p.Results.PixelSize);
            end
            if ~isnan(p.Results.FrameInterval)
                tw.setFrameInterval(p.Results.FrameInterval);
            end

            tw.readData(csvFile);
            try
                tw.loadNucleusMask();
                tw.cull_tracks();
                tw.remove_relative_motion();
            catch
                warning('no nucleus file or cullable tracks');
            end

            obj = obj.addWrapper(tw, 'FilePath', csvFile, ...
                'FileID', fileID, 'Condition', p.Results.Condition);
        end

        %% Add from a processed SMD object
        function obj = addFromSMD(obj, smd, varargin)
            % ADDFROMSMD  Create a TrajectoryWrapper from a processed SMD object and add.
            %
            %   tc.addFromSMD(smd)
            %   tc.addFromSMD(smd, 'FileID','run1', 'Condition','Control')
            %
            %   The SMD object must already be processed (localize → track →
            %   get_roi → cull_tracks). No pipeline steps are run here.

            p = inputParser;
            addParameter(p, 'FileID',    '');
            addParameter(p, 'Condition', '');
            parse(p, varargin{:});

            fileID = p.Results.FileID;
            if isempty(fileID)
                [~, fileID, ~] = fileparts(smd.fname);
            end

            tw = TrajectoryWrapper.fromSMD(smd);

            obj = obj.addWrapper(tw, 'FilePath', smd.fname, ...
                'FileID', fileID, 'Condition', p.Results.Condition);
        end

        %% Add existing TrajectoryWrapper
        function obj = addWrapper(obj, tw, varargin)
            if ~isa(tw, 'TrajectoryWrapper')
                error('Input must be a TrajectoryWrapper instance.');
            end
            p = inputParser;
            addParameter(p, 'FileID',   num2str(length(obj.Wrappers) + 1));
            addParameter(p, 'FilePath', '');
            addParameter(p, 'Condition','');
            parse(p, varargin{:});

            obj.Wrappers{end+1} = tw;
            obj.Metadata = [obj.Metadata; {p.Results.FileID, p.Results.FilePath, p.Results.Condition}];
            obj.invalidateCollections();
        end

        %% Batch load from file lists
        function obj = loadFromFiles(obj, csvFiles, maskFiles, checkFunc, conditions)
            if length(csvFiles) ~= length(maskFiles)
                error('CSV and mask file lists must match in length.');
            end
            if nargin < 5 || isempty(conditions)
                conditions = repmat({''}, length(csvFiles), 1);
            elseif ischar(conditions) || isstring(conditions)
                conditions = repmat({conditions}, length(csvFiles), 1);
            end
            for i = 1:length(csvFiles)
                obj.addFromFile(csvFiles{i}, maskFiles{i}, checkFunc, ...
                    'FileID', num2str(i), 'Condition', conditions{i});
            end
        end

        %% Read experiment params from TOML
        function obj = ReadParams(obj, fileID)
            % READPARAMS  Load TOML parameters into Params and propagate to wrappers.
            %
            %   Recognized top-level TOML keys (propagated to wrappers):
            %     pixel_size      → sets PixelSize on CSV-sourced wrappers
            %     frame_interval  → sets FrameInterval on wrappers where it is NaN
            %   All keys are also stored in obj.Params for analysis code to read.

            config     = toml.read(fileID);
            pars       = toml.map_to_struct(config);
            obj.Params = pars;

            % Propagate pixel_size / frame_interval to existing wrappers
            for w = 1:length(obj.Wrappers)
                tw = obj.Wrappers{w};
                if isfield(pars, 'pixel_size') && strcmp(tw.getSourceType(), 'csv')
                    % Only override for CSV wrappers; SMD wrappers derive px from camera cal
                    tw.setPixelSize(pars.pixel_size);
                end
                if isfield(pars, 'frame_interval') && isnan(tw.getFrameInterval())
                    tw.setFrameInterval(pars.frame_interval);
                end
            end
            obj.invalidateCollections();
        end

        %% Set frame interval on all wrappers
        function obj = setFrameIntervalAll(obj, dt)
            % SETFRAMEINTERVALALL  Set FrameInterval (s) on every wrapper in the collection.
            %
            %   tc.setFrameIntervalAll(0.02)
            %
            %   Overrides existing values on all wrappers (including SMD-sourced ones).
            %   Call this after all addFromFile / addFromSMD calls are complete.
            for w = 1:length(obj.Wrappers)
                obj.Wrappers{w}.setFrameInterval(dt);
            end
            obj.invalidateCollections();
        end

        %% Lazy accessors
        function tracks = getAllRawTracks(obj)
            if ~obj.IsAllRawCollected, obj.collectAllTracks(); end
            tracks = obj.AllRawTracks;
        end

        function tracks = getAllCulledTracks(obj)
            if ~obj.IsAllCulledCollected, obj.collectAllTracks(); end
            tracks = obj.AllCulledTracks;
        end

        function tracks = getAllRelativeTracks(obj)
            if ~obj.IsAllRelativeCollected, obj.collectAllTracks(); end
            tracks = obj.AllRelativeTracks;
        end

        %% RL decomposition (Richardson-Lucy MSD distribution)
        function results = getRLDecomposition(obj, varargin)
            % GETRLDECOMPOSITION  Richardson-Lucy decomposition of the MSD distribution.
            %
            %   results = tc.getRLDecomposition('LagTime', 4, 'String', 'H2B')
            %   results = tc.getRLDecomposition('TrackType', 'relative', ...)
            %
            %   Runs RL_HILO on the selected track population and stores results in
            %   tc.RLResults.  Re-call with the same arguments returns cached results;
            %   use 'ForceRecompute' to re-run (also required when switching TrackType).
            %
            %   Options:
            %     'TrackType'        - 'culled' (default) | 'relative' | 'raw'
            %                          Use 'relative' to analyse motion after cell-
            %                          motion subtraction and compare state fractions
            %                          and diffusion parameters with the culled result.
            %     'LagTime'          - frame lag for van Hove computation (default 4)
            %     'String'           - label for figure titles (default 'Not specified')
            %     'Condition'        - restrict to one condition (default: all)
            %     'ExposureFraction' - Te/dt for fBM model inside RL (default 1)
            %     'MinStepVar'       - minimum mean squared step (µm²) to include a
            %                          track; drops stuck particles before RL runs.
            %                          Set to ~4*sigma^2 (default: 0 = no filter).
            %     'MinTrackLength'   - minimum track length in frames (default: 7)
            %     'MinGroupSize'     - minimum tracks per state to run fBM fit (default: 50)
            %     'SubtrackLength'   - subtrack length for per-state fBM MLE (default: 20)
            %     'CIMethod'         - CI method for per-state fBM MLE (default: 'none')
            %     'ForceRecompute'   - logical; default false

            varargin = unpack_opts(varargin{:});
            p = inputParser;
            addParameter(p, 'TrackType',        'culled',        @(x) ismember(lower(char(x)), {'culled','raw','relative'}));
            addParameter(p, 'LagTime',          4,               @isnumeric);
            addParameter(p, 'String',            'Not specified', @(x) ischar(x)||isstring(x));
            addParameter(p, 'Condition',         [],              @(x) ischar(x)||isstring(x)||iscell(x));
            addParameter(p, 'ExposureFraction',  1,               @isnumeric);
            addParameter(p, 'MinStepVar',        0,               @isnumeric);
            addParameter(p, 'MinTrackLength',    7,               @isnumeric);
            addParameter(p, 'MinGroupSize',      50,              @isnumeric);
            addParameter(p, 'SubtrackLength',    20,              @isnumeric);
            addParameter(p, 'CIMethod',          'none',          @ischar);
            addParameter(p, 'ForceRecompute',    false,           @islogical);
            parse(p, varargin{:});
            o = p.Results;
            trackType = lower(char(o.TrackType));

            % Cache is valid only if the same TrackType was used previously
            cachedType = '';
            if obj.IsRLComputed && isfield(obj.RLResults, 'track_type')
                cachedType = obj.RLResults.track_type;
            end
            if obj.IsRLComputed && ~o.ForceRecompute && strcmp(cachedType, trackType)
                results = obj.RLResults;
                return;
            end
            if obj.IsRLComputed && ~o.ForceRecompute && ~strcmp(cachedType, trackType)
                warning('getRLDecomposition: cached result used TrackType=''%s'' but ''%s'' requested. Pass ''ForceRecompute'',true to rerun.', ...
                    cachedType, trackType);
                results = obj.RLResults;
                return;
            end

            % Retrieve tracks for the requested type
            dt = obj.getFrameIntervalForCondition(o.Condition);
            switch trackType
                case 'culled'
                    if isempty(o.Condition)
                        tracks = obj.getAllCulledTracks();
                    else
                        tracks = obj.getTracksByCondition(o.Condition);
                    end
                case 'relative'
                    if isempty(o.Condition)
                        tracks = obj.getAllRelativeTracks();
                    else
                        tracks = obj.getRelativeTracksByCondition(o.Condition);
                    end
                    if isempty(tracks)
                        error('getRLDecomposition: no relative tracks found. Run remove_relative_motion() on your wrappers first.');
                    end
                case 'raw'
                    if isempty(o.Condition)
                        tracks = obj.getAllRawTracks();
                    else
                        tracks = obj.getTracksByCondition(o.Condition, 'UseRaw', true);
                    end
            end

            if isempty(tracks)
                error('getRLDecomposition: no tracks available.');
            end

            % Length filter
            Lmin     = round(o.MinTrackLength);
            tracklen = cellfun('size', tracks, 1);
            keep     = tracklen >= Lmin;

            % Stuck-particle pre-filter
            if o.MinStepVar > 0
                mobile = cellfun(@(t) mean(diff(t(:,1)).^2 + diff(t(:,2)).^2) >= o.MinStepVar, tracks);
                keep   = keep & mobile;
                fprintf('getRLDecomposition: %d/%d tracks removed by MinStepVar filter (< %.4g µm²)\n', ...
                    sum(~mobile & tracklen >= Lmin), sum(tracklen >= Lmin), o.MinStepVar);
            end

            % keep_idx maps position-in-filtered-array → position-in-full-array.
            % classified_tracks indices from RL_HILO are into the filtered array,
            % so getFBMByRLState translates them through this map.
            keep_idx = find(keep);

            X = cell(1, numel(tracks));
            for itrack = 1:numel(tracks)
                X{itrack} = tracks{itrack}(:, 1:2);
            end
            tmp{1} = X(keep);

            fprintf('getRLDecomposition: %d tracks entering RL_HILO\n', numel(keep_idx));

            [M, ~, computed_quantities] = RL_HILO(tmp, o.String, o.LagTime, ...
                'dt',               dt, ...
                'ExposureFraction', o.ExposureFraction, ...
                'MinGroupSize',     o.MinGroupSize, ...
                'SubtrackLength',   o.SubtrackLength, ...
                'CIMethod',         o.CIMethod);

            obj.RLResults.M                   = M;
            obj.RLResults.computed_quantities = computed_quantities;
            obj.RLResults.dt                  = dt;
            obj.RLResults.frac                = o.ExposureFraction;
            obj.RLResults.track_type          = trackType;
            obj.RLResults.track_index_map     = keep_idx;   % filtered→full index translation
            obj.IsRLComputed = true;
            results = obj.RLResults;
        end

        %% Ensemble MSD with bootstrap fitting
        function results = getMSD(obj, varargin)
            % GETMSD  Ensemble-averaged MSD with bootstrap power-law fits.
            %
            %   results = tc.getMSD()
            %   results = tc.getMSD('TrackType',    'culled')   % 'culled'|'relative'|'raw'
            %   results = tc.getMSD('MaxLag',       25)         % max lag in frames
            %   results = tc.getMSD('NumBootstrap', 50)         % bootstrap samples
            %   results = tc.getMSD('ExposureTime', 0.01)       % camera exposure (s)
            %   results = tc.getMSD('Condition',    'Control')
            %   results = tc.getMSD('ForceRecompute', true)
            %
            %   ExposureTime defaults to Params.exposure_time_s if set via TOML,
            %   otherwise to the frame interval (frac = 1, i.e. stroboscopic).
            %
            %   Results struct fields:
            %     .ensemble_mean  - [MaxLag x 1] ensemble-averaged MSD (µm²)
            %     .lag_axis       - [MaxLag x 1] lag times in seconds
            %     .lag_frames     - [MaxLag x 1] lag indices (1 … MaxLag)
            %     .boot_means     - [MaxLag x NumBootstrap] per-bootstrap means
            %     .boot_fits      - [NumBootstrap x 3] fit params [G, sigma_sq, alpha]
            %     .dt             - frame interval (s)
            %     .frac           - exposure / frame_interval ratio
            %     .track_type     - which track set was used
            %
            %   Note: the fBM fit is performed over the first floor(MaxLag/2) lag
            %   points only. The fit curve in plotMSD is extrapolated to the full
            %   lag range for display.

            varargin = unpack_opts(varargin{:});
            p = inputParser;
            addParameter(p, 'TrackType',      'culled', @(x) ischar(x)||isstring(x));
            addParameter(p, 'MaxLag',         25,       @isnumeric);
            addParameter(p, 'NumBootstrap',   50,       @isnumeric);
            addParameter(p, 'ExposureTime',   NaN,      @isnumeric);
            addParameter(p, 'MinLength',       1,       @isnumeric);
            addParameter(p, 'MinMeanSqStep',  0,       @isnumeric);
            addParameter(p, 'Condition',      [],       @(x) ischar(x)||isstring(x)||iscell(x));
            addParameter(p, 'ForceRecompute', false,    @islogical);
            parse(p, varargin{:});
            o = p.Results;

            if p.Results.ForceRecompute
                obj.IsMSDComputed = false;
                obj.MSDResults    = struct();
            end

            if ~obj.IsMSDComputed

                % --- Get tracks -------------------------------------------
                trackType = lower(char(o.TrackType));
                switch trackType
                    case 'culled'
                        if isempty(o.Condition)
                            tracks = obj.getAllCulledTracks();
                        else
                            tracks = obj.getTracksByCondition(o.Condition);
                        end
                    case 'relative'
                        if isempty(o.Condition)
                            tracks = obj.getAllRelativeTracks();
                        else
                            tracks = obj.getRelativeTracksByCondition(o.Condition);
                        end
                    case 'raw'
                        if isempty(o.Condition)
                            tracks = obj.getAllRawTracks();
                        else
                            tracks = obj.getTracksByCondition(o.Condition, 'UseRaw', true);
                        end
                    otherwise
                        error('getMSD: TrackType must be ''culled'', ''relative'', or ''raw''.');
                end

                if isempty(tracks)
                    error('getMSD: no tracks available (TrackType=''%s'').', trackType);
                end

                % --- MinLength filter -------------------------------------
                if o.MinLength > 1
                    tracks = tracks(cellfun(@(t) size(t,1), tracks) >= o.MinLength);
                    if isempty(tracks)
                        error('getMSD: no tracks of length >= %d (TrackType=''%s'').', ...
                              o.MinLength, trackType);
                    end
                end

                % --- Immobility filter ------------------------------------
                if o.MinMeanSqStep > 0
                    n_before = numel(tracks);
                    tracks   = filter_mobile(tracks, o.MinMeanSqStep);
                    if isempty(tracks)
                        error('getMSD: no mobile tracks (MinMeanSqStep=%.4g).', o.MinMeanSqStep);
                    end
                    fprintf('getMSD: removed %d immobile tracks (MinMeanSqStep=%.4g µm²), %d remain.\n', ...
                            n_before - numel(tracks), o.MinMeanSqStep, numel(tracks));
                end

                % --- Physical parameters ----------------------------------
                dt = obj.getFrameIntervalForCondition(o.Condition);

                expT = o.ExposureTime;
                if isnan(expT) && isfield(obj.Params, 'exposure_time_s')
                    expT = obj.Params.exposure_time_s;
                end
                if isnan(expT)
                    frac = 1;   % assume stroboscopic (exposure = frame interval)
                    warning('getMSD:noExposureTime', ...
                        ['ExposureTime not set; assuming frac=1 (stroboscopic). ' ...
                         'Pass ExposureTime or set exposure_time_s in your TOML.']);
                else
                    frac = expT / dt;
                end

                MaxLag = o.MaxLag;

                % --- Gap-fill and per-track MSD ---------------------------
                nTracks      = numel(tracks);
                ensemble_msd  = cell(nTracks, 1);
                ensemble_nlags = cell(nTracks, 1);

                for ii = 1:nTracks
                    t_raw = tracks{ii};
                    x = t_raw(:, 1);
                    y = t_raw(:, 2);
                    if size(t_raw, 2) >= 3
                        t = t_raw(:, 3) - t_raw(1, 3) + 1;   % 1-based frame index
                    else
                        t = (1:length(x))';
                    end
                    [t_u, x_u] = fillGapsWithNaN(t, x);
                    [~,   y_u] = fillGapsWithNaN(t, y);
                    [msd_i, nlags_i] = simple_msd([x_u, y_u], t_u);
                    ensemble_msd{ii}   = msd_i;
                    ensemble_nlags{ii} = nlags_i;
                end

                % --- Ensemble mean ----------------------------------------
                msdlen = cellfun(@length, ensemble_msd);
                N      = max(msdlen);
                e_sum  = zeros(N, 1);   % weighted sum of MSDs
                n_sum  = zeros(N, 1);   % total observation count at each lag

                for ii = 1:nTracks
                    L = msdlen(ii);
                    e_sum(1:L) = e_sum(1:L) + ensemble_msd{ii};
                    n_sum(1:L) = n_sum(1:L) + ensemble_nlags{ii};
                end
                n_sum(n_sum == 0) = NaN;
                ensemble_mean_full = e_sum ./ n_sum;

                % Trim to MaxLag
                nLags         = min(MaxLag, N);
                ensemble_mean = ensemble_mean_full(1:nLags);
                lag_frames    = (1:nLags)';
                lag_axis      = lag_frames * dt;

                % --- Bootstrap -------------------------------------------
                nBoot      = o.NumBootstrap;
                boot_means = zeros(nLags, nBoot);
                boot_fits  = zeros(nBoot, 3);

                lsqOpts = optimoptions('lsqnonlin', ...
                    'MaxFunctionEvaluations', 10000, ...
                    'Display', 'none', ...
                    'Algorithm', 'levenberg-marquardt');

                x0  = [1e-4, 4e-4, 0.5];
                lb  = [1e-5, 1e-6, 0.1];
                ub  = [5,    4e-3, 1.5];

                for isamp = 1:nBoot
                    idx  = datasample(1:nTracks, nTracks);
                    bs_sum  = zeros(N, 1);
                    bs_n    = zeros(N, 1);
                    for jj = 1:nTracks
                        ii = idx(jj);
                        L  = msdlen(ii);
                        bs_sum(1:L) = bs_sum(1:L) + ensemble_msd{ii};
                        bs_n(1:L)   = bs_n(1:L)   + ensemble_nlags{ii};
                    end
                    bs_n(bs_n == 0) = NaN;
                    bm = bs_sum ./ bs_n;
                    boot_means(:, isamp) = bm(1:nLags);

                    fitLags = 1:max(1, floor(nLags/2));
                    bm_fit = bm(fitLags);
                    valid  = ~isnan(bm_fit) & bm_fit > 0;
                    if sum(valid) >= 3
                        try
                            fp = lsqnonlin( ...
                                @(x) my_fun(x, lag_frames(fitLags(valid)), bm_fit(valid), dt, frac), ...
                                x0, lb, ub, lsqOpts);
                            boot_fits(isamp, :) = fp;
                        catch
                            boot_fits(isamp, :) = NaN;
                        end
                    else
                        boot_fits(isamp, :) = NaN;
                    end
                end

                % --- Store ------------------------------------------------
                obj.MSDResults.ensemble_mean = ensemble_mean;
                obj.MSDResults.lag_axis      = lag_axis;
                obj.MSDResults.lag_frames    = lag_frames;
                obj.MSDResults.boot_means    = boot_means;
                obj.MSDResults.boot_fits     = boot_fits;
                obj.MSDResults.dt            = dt;
                obj.MSDResults.frac          = frac;
                obj.MSDResults.track_type    = trackType;
                obj.MSDResults.n_tracks      = nTracks;
                obj.MSDResults.min_length    = o.MinLength;
                obj.IsMSDComputed = true;

                fprintf('getMSD: %d tracks (MinLength=%d), dt=%.4fs, frac=%.3f, MaxLag=%d, %d bootstrap samples.\n', ...
                    nTracks, o.MinLength, dt, frac, nLags, nBoot);
            end
            results = obj.MSDResults;
        end

        %% Bayesian diffusivity
        function results = getBayesianDiffusivity(obj, varargin)
            % GETBAYESIANDIFFUSIVITY  Bayesian (pEMv2) diffusivity estimation.
            %
            %   results = tc.getBayesianDiffusivity()
            %   results = tc.getBayesianDiffusivity('Condition','Control')
            %
            % FrameInterval is read per-condition from wrappers.  If any
            % wrapper has FrameInterval=NaN this will error with a clear message.

            varargin = unpack_opts(varargin{:});
            p = inputParser;
            addParameter(p, 'Condition',      [], @(x) ischar(x)||isstring(x)||iscell(x));
            addParameter(p, 'ForceRecompute', [], @islogical);
            addParameter(p, 'MinLength',       6, @isnumeric);
            addParameter(p, 'SplitLength',     7, @isnumeric);
            addParameter(p, 'NumFeatures',    [], @isnumeric);  % default: splitLength-1, max 10
            parse(p, varargin{:});
            opts = p.Results;

            isSubset = ~isempty(opts.Condition);

            if isempty(opts.ForceRecompute) && obj.IspEMComputed
                force = false;
            elseif isempty(opts.ForceRecompute) && ~obj.IspEMComputed
                force = true;
            elseif ~isempty(opts.ForceRecompute) && opts.ForceRecompute
                force = true;
            end
            if isSubset
                force = true;
            end

            fprintf('\n=== Bayesian Diffusivity Analysis ===\n');
            fprintf('Timestamp: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

            if isSubset
                condStr = strjoin(cellstr(opts.Condition), ', ');
                fprintf('Restricted to condition(s): %s\n', condStr);
            else
                fprintf('Analyzing full collection (%d FOVs)\n', length(obj.Wrappers));
            end

            if isSubset
                tracks = obj.getTracksByCondition(opts.Condition);
                dt     = obj.getFrameIntervalForCondition(opts.Condition);
            else
                tracks = obj.getAllCulledTracks();
                dt     = obj.getFrameIntervalForCondition([]);
            end

            if isempty(tracks)
                warning('No tracks available for Bayesian analysis.');
                results = struct();
                fprintf('=====================================\n\n');
                return;
            end

            if force
                fprintf('Analyzing (%d tracks)\n', numel(tracks));
                splitLength = opts.SplitLength;
                if isempty(opts.NumFeatures)
                    numFeatures = min(splitLength - 1, 10);
                else
                    numFeatures = opts.NumFeatures;
                end
                dE             = 0.01;
                minStates      = 1;
                maxStates      = 8;
                numReinitialize = 20;
                numPerturb      = 200;
                maxiter         = 10000;
                convergence     = 1e-7;
                lambda          = 0.0;

                k       = 1;
                X       = {};
                trackID = [];
                for itrack = 1:numel(tracks)
                    if size(tracks{itrack}, 1) > opts.MinLength
                        X{k}       = tracks{itrack}(:, 1:2);   % already in µm
                        trackID(k) = k;
                        k = k + 1;
                    end
                end

                [splitX, splitIndex] = SplitTracks(X, splitLength);

                trackInfo.numberOfTracks = length(splitX);
                trackInfo.dimensions     = size(splitX{1}, 2);
                trackInfo.numFeatures    = numFeatures;
                trackInfo.splitLength    = splitLength;
                trackInfo.splitIndex     = splitIndex;
                trackInfo.dt             = dt;
                trackInfo.R              = 1/6 * dE / dt;
                trackInfo.lambda         = lambda;

                params.minStates       = minStates;
                params.maxStates       = maxStates;
                params.numFeatures     = numFeatures;
                params.numReinitialize = numReinitialize;
                params.numPerturbation = numPerturb;
                params.converged       = convergence;
                params.maxiter         = maxiter;
                params.verbose         = 1;

                deltaX = cell(trackInfo.numberOfTracks, 1);
                for i = 1:trackInfo.numberOfTracks
                    deltaX{i} = diff(splitX{i});
                end

                [trackInfo.vacf_exp, trackInfo.xbar_exp] = ...
                    CovarianceProperties(deltaX, numFeatures);

                tic;
                results = pEMv2_SPT(deltaX, trackInfo, params);
                toc;

                % Store inputs needed for bootstrap CI (computed separately via
                % computeBootstrapCI — not run here so a parallel pool failure
                % cannot prevent pEM results from being saved).
                if ~isSubset
                    obj.pEMBootstrapInputs.deltaX     = deltaX;
                    obj.pEMBootstrapInputs.rawResults = results;
                    obj.pEMBootstrapInputs.trackInfo  = trackInfo;
                    obj.pEMBootstrapInputs.params     = params;
                    obj.IsBootstrapComputed           = false;
                end

                vacf              = results.optimalVacf;
                results.optimalD  = (vacf(:,1) + 2*vacf(:,2)) / (2*trackInfo.dt);
                results.optimalS  = sqrt((vacf(:,1) - 2*results.optimalD*dt*(1-2*trackInfo.R)) / 2);

                [D, idx] = sort(results.optimalD);

                pEMTable = table;
                pEMTable.numRawTracks{1}   = numel(X);
                pEMTable.numSplitTracks(1) = results.trackInfo.numberOfTracks;
                pEMTable.trackInfo{1}      = results.trackInfo;
                pEMTable.params{1}         = results.params;
                pEMTable.optimalSize(1)    = results.optimalSize(1);
                pEMTable.optimalD{1}       = D;
                pEMTable.optimalS{1}       = results.optimalS(idx);
                pEMTable.optimalVacf{1}    = results.optimalVacf(idx, :);
                pEMTable.optimalP{1}       = results.optimalP(idx);
                pEMTable.optimalL(1)       = results.optimalL;
                pEMTable.BIC{1}            = results.BIC(idx);
                pEMTable.splitX{1}         = splitX;
                pEMTable.trackID{1}        = trackID;
                pEMTable.splitID{1}        = results.trackInfo.splitIndex;
                pEMTable.posteriorProb{1}  = results.posteriorProb(:, idx);

                posteriorProb              = pEMTable.posteriorProb{1};
                [maxPosteriorProb, optimalState] = max(posteriorProb, [], 2);

                vacf    = pEMTable.trackInfo{1}.vacf_exp;
                dt_tbl  = pEMTable.trackInfo{1}.dt;
                track_D = mean(vacf(:,1,:) + 2*vacf(:,2,:), 3) / 2 / dt_tbl;

                pEMTable.maxPosteriorProb{1} = maxPosteriorProb;
                pEMTable.optimalState{1}     = optimalState;
                pEMTable.track_D{1}          = track_D;

                maxPosteriorProb = pEMTable.maxPosteriorProb{1};
                optimalState     = pEMTable.optimalState{1};

                proportions = table;
                for i = 1:pEMTable.optimalSize
                    proportions.State(i)      = i;
                    proportions.Proportion(i) = length(find(optimalState==i)) / length(optimalState);
                    proportions.OptimalD(i)   = pEMTable.optimalD{1}(i);
                    proportions.OptimalS(i)   = pEMTable.optimalS{1}(i);
                end

                fprintf('Population fraction (global):\n');
                proportions

                if ~isSubset
                    obj.pEMParams.trackInfo    = trackInfo;
                    obj.pEMParams.params       = params;
                    obj.IspEMComputed          = true;
                    obj.pEMResults.pEMTable    = pEMTable;
                    obj.pEMResults.proportions = proportions;
                end
            end
            results = obj.pEMResults;
        end

        %% Population fBM MLE  (K, alpha, sigma)
        function results = getFBMParameters(obj, varargin)
            % GETFBMPARAMETERS  Population MLE of fBM parameters via Toeplitz likelihood.
            %
            %   results = tc.getFBMParameters()
            %   results = tc.getFBMParameters('TrackType','culled', 'CIMethod','profile', ...)
            %
            %   Fits a single (K, alpha, sigma) to all tracks by maximising the
            %   exact multivariate-Gaussian log-likelihood of the displacement
            %   time series (Backlund et al. Phys. Rev. E 2015 covariance model).
            %   The sufficient-statistic formulation makes the per-evaluation cost
            %   O(L^3) regardless of track count, so 100 K tracks cost no more than 1 K.
            %
            %   Each track longer than SubtrackLength is split into non-overlapping
            %   segments of exactly SubtrackLength frames; the final short remainder
            %   is kept if it is >= MinSubtrackLength, discarded otherwise.
            %
            %   Options (passed through to fitFBM_MLE):
            %     'TrackType'         - 'culled' (default) | 'raw' | 'relative'
            %                           NOTE: 'relative' removes global cell motion
            %                           but introduces a bias if residual motion
            %                           correlations remain; use with caution.
            %     'Condition'         - filter by condition string
            %     'ExposureFraction'  - Te/dt  (default: 1 = stroboscopic)
            %     'SubtrackLength'    - frames per subtrack (default: 20)
            %     'MinSubtrackLength' - minimum subtrack length kept (default: 10)
            %     'CIMethod'          - 'profile'|'profile+mcmc'|'mcmc'|'hessian'|'none'
            %                           (default: 'profile')
            %     'ProfilePoints'     - grid points per profile (default: 60)
            %     'MCMCSamples'       - posterior samples if using MCMC (default: 5000)
            %     'ForceRecompute'    - rerun even if cached (default: false)
            %
            %   Result fields: K, Ka (=K*alpha), alpha, sigma, loglik,
            %                  n_subtracks, n_displacements, CI (if requested).

            varargin = unpack_opts(varargin{:});
            p = inputParser;
            addParameter(p, 'TrackType',         'culled', @(x) ismember(lower(char(x)), {'culled','raw','relative'}));
            addParameter(p, 'Condition',         [],       @(x) ischar(x)||isstring(x)||iscell(x)||isempty(x));
            addParameter(p, 'ExposureFraction',  1,        @isnumeric);
            addParameter(p, 'SubtrackLength',    20,       @isnumeric);
            addParameter(p, 'MinSubtrackLength', 10,       @isnumeric);
            addParameter(p, 'MinMeanSqStep',     0,        @isnumeric);
            addParameter(p, 'CIMethod',          'profile',@ischar);
            addParameter(p, 'ProfilePoints',     60,       @isnumeric);
            addParameter(p, 'MCMCSamples',       5000,     @isnumeric);
            addParameter(p, 'ForceRecompute',    false,    @islogical);
            parse(p, varargin{:});
            o = p.Results;

            trackType = lower(char(o.TrackType));
            isSubset  = ~isempty(o.Condition) || ~strcmp(trackType, 'culled');

            % Return cached result only if all fitting options match what was previously used.
            cached_ci_method = '';
            if obj.IsFBMComputed && isfield(obj.FBMResults, 'CI') && isstruct(obj.FBMResults.CI)
                cached_ci_method = obj.FBMResults.CI.method;
            elseif obj.IsFBMComputed && (~isfield(obj.FBMResults, 'CI') || isempty(obj.FBMResults.CI))
                cached_ci_method = 'none';
            end
            ci_match   = strcmpi(cached_ci_method, o.CIMethod);
            frac_match = obj.IsFBMComputed && isfield(obj.FBMResults, 'frac') && ...
                         obj.FBMResults.frac == o.ExposureFraction;
            lsub_match = obj.IsFBMComputed && isfield(obj.FBMResults, 'subtrack_length') && ...
                         obj.FBMResults.subtrack_length == o.SubtrackLength;
            lmin_match = obj.IsFBMComputed && isfield(obj.FBMResults, 'min_subtrack_length') && ...
                         obj.FBMResults.min_subtrack_length == o.MinSubtrackLength;
            mss_match  = obj.IsFBMComputed && isfield(obj.FBMResults, 'min_mean_sq_step') && ...
                         obj.FBMResults.min_mean_sq_step == o.MinMeanSqStep;

            if ~o.ForceRecompute && obj.IsFBMComputed && ~isSubset && ...
               ci_match && frac_match && lsub_match && lmin_match && mss_match
                results = obj.FBMResults;
                return;
            end

            [tracks, dt] = obj.getTracksForFBM(trackType, o.Condition);

            if isempty(tracks)
                warning('getFBMParameters: no tracks found (TrackType=''%s'').', trackType);
                results = struct();
                return;
            end

            % --- Immobility filter ----------------------------------------
            if o.MinMeanSqStep > 0
                n_before = numel(tracks);
                tracks   = filter_mobile(tracks, o.MinMeanSqStep);
                if isempty(tracks)
                    error('getFBMParameters: no mobile tracks (MinMeanSqStep=%.4g).', o.MinMeanSqStep);
                end
                fprintf('getFBMParameters: removed %d immobile tracks (MinMeanSqStep=%.4g µm²), %d remain.\n', ...
                        n_before - numel(tracks), o.MinMeanSqStep, numel(tracks));
            end

            fprintf('getFBMParameters: %d %s tracks, dt=%.4f s, CIMethod=%s\n', ...
                numel(tracks), trackType, dt, o.CIMethod);

            [fp, ~, ci] = fitFBM_MLE(tracks, dt, ...
                'ExposureFraction',  o.ExposureFraction, ...
                'SubtrackLength',    o.SubtrackLength, ...
                'MinSubtrackLength', o.MinSubtrackLength, ...
                'CIMethod',          o.CIMethod, ...
                'ProfilePoints',     o.ProfilePoints, ...
                'MCMCSamples',       o.MCMCSamples);
            results                    = fp;
            results.CI                 = ci;
            results.track_type         = trackType;
            results.frac               = o.ExposureFraction;
            results.subtrack_length    = o.SubtrackLength;
            results.min_subtrack_length = o.MinSubtrackLength;
            results.min_mean_sq_step   = o.MinMeanSqStep;

            fprintf('  K=%.4g µm²/s^a  alpha=%.3f  sigma=%.3f µm  loglik=%.1f  (%d subtracks)\n', ...
                results.K, results.alpha, results.sigma, results.loglik, results.n_subtracks);

            if ~isSubset
                obj.FBMResults    = results;
                obj.IsFBMComputed = true;
            end
        end

        %% Per-track (K, alpha) distribution
        function results = getFBMAlphaDistribution(obj, varargin)
            % GETFBMALPHADISTRIBUTION  Per-track fBM parameter estimates.
            %
            %   results = tc.getFBMAlphaDistribution()
            %   results = tc.getFBMAlphaDistribution('TrackType','culled', 'FixSigma',true)
            %
            %   Fits (K, alpha) independently to each track.  With 'FixSigma' true
            %   (default), sigma is held fixed at the population estimate from
            %   getFBMParameters(), reducing the per-track problem from 3D to 2D
            %   and substantially reducing noise in the alpha estimates.
            %
            %   With only ~38 scalar observations per 20-frame track, per-track
            %   alpha estimates have inherent variance (~0.2-0.4 s.d.); the
            %   observed distribution is a convolution of the true biological
            %   distribution with this fitting noise.
            %
            %   Options:
            %     'TrackType'        - 'culled' (default) | 'raw' | 'relative'
            %                          NOTE: 'relative' tracks have cell motion
            %                          subtracted; the fBM covariance assumes
            %                          stationary increments, so any residual
            %                          correlated motion biases alpha upward.
            %                          Use 'culled' unless you have verified the
            %                          subtraction is clean.
            %     'Condition'        - filter by condition
            %     'ExposureFraction' - Te/dt (default: 1)
            %     'FixSigma'         - fix sigma to population value (default: true)
            %     'SigmaValue'       - override sigma (µm); bypasses auto-estimation
            %     'MinTrackLength'     - minimum frames per track (default: 15)
            %     'MaxSubtrackLength'  - cap the Toeplitz matrix size by splitting
            %                            long tracks into non-overlapping subtracks
            %                            of this length and pooling their sufficient
            %                            statistics.  Set to the same value as
            %                            SubtrackLength in getFBMParameters (e.g. 20)
            %                            so per-track and population fits use the
            %                            same lag range.  Default: 0 = full track.
            %     'MinStepVar'         - minimum mean squared step size (µm²) to
            %                            include a track (default: 0 = no filter).
            %                            Use e.g. 4*sigma^2 to drop stuck particles
            %                            before fitting; these tracks collapse alpha
            %                            toward zero and bias DACF diagnostics.
            %     'MinAlpha'           - post-fit minimum alpha threshold; tracks
            %                            with fitted alpha < MinAlpha are removed
            %                            from the returned results (default: 0 = keep
            %                            all).  Use e.g. 0.1 to discard residual
            %                            near-immobile particles.
            %     'ForceRecompute'     - rerun even if cached (default: false)
            %     'Verbose'            - print progress (default: false)
            %
            %   Result fields:
            %     .alpha        Nx1 per-track anomalous exponent
            %     .K            Nx1 per-track K (µm²/s^alpha)
            %     .Ka           Nx1 K*alpha
            %     .track_length Nx1 frames in each track
            %     .sigma_fixed  scalar sigma used (NaN if fitted jointly)
            %     .track_type   string
            %     .n_removed    number of tracks removed by MinAlpha filter

            varargin = unpack_opts(varargin{:});
            p = inputParser;
            addParameter(p, 'TrackType',        'culled', @(x) ismember(lower(char(x)), {'culled','raw','relative'}));
            addParameter(p, 'Condition',        [],    @(x) ischar(x)||isstring(x)||iscell(x)||isempty(x));
            addParameter(p, 'ExposureFraction',  1,    @isnumeric);
            addParameter(p, 'FixSigma',       true,    @islogical);
            addParameter(p, 'SigmaValue',      NaN,    @isnumeric);
            addParameter(p, 'MinTrackLength',    15,    @isnumeric);
            addParameter(p, 'MaxSubtrackLength',  0,    @isnumeric);
            addParameter(p, 'MinStepVar',          0,    @isnumeric);
            addParameter(p, 'MinAlpha',            0,    @isnumeric);
            addParameter(p, 'ForceRecompute',  false,   @islogical);
            addParameter(p, 'Verbose',         false,   @islogical);
            parse(p, varargin{:});
            o = p.Results;

            trackType = lower(char(o.TrackType));
            isSubset  = ~isempty(o.Condition) || ~strcmp(trackType, 'culled');

            if ~o.ForceRecompute && obj.IsFBMAlphaComputed && ~isSubset
                results = obj.FBMAlphaResults;
                return;
            end

            [tracks, dt] = obj.getTracksForFBM(trackType, o.Condition);

            if isempty(tracks)
                warning('getFBMAlphaDistribution: no tracks found (TrackType=''%s'').', trackType);
                results = struct();
                return;
            end

            % Determine sigma to fix
            sigma_fixed = NaN;
            if o.FixSigma
                if ~isnan(o.SigmaValue)
                    sigma_fixed = o.SigmaValue;
                elseif obj.IsFBMComputed && isfield(obj.FBMResults, 'sigma') && ...
                        strcmp(obj.FBMResults.track_type, trackType)
                    sigma_fixed = obj.FBMResults.sigma;
                else
                    fprintf('getFBMAlphaDistribution: running population MLE to get sigma...\n');
                    pop = obj.getFBMParameters('TrackType', trackType, ...
                        'ExposureFraction', o.ExposureFraction, 'CIMethod', 'none');
                    sigma_fixed = pop.sigma;
                end
                fprintf('getFBMAlphaDistribution: sigma fixed at %.4f µm\n', sigma_fixed);
            end

            results = fitFBM_pertracks(tracks, dt, ...
                'ExposureFraction',   o.ExposureFraction, ...
                'MinTrackLength',     o.MinTrackLength, ...
                'MaxSubtrackLength',  o.MaxSubtrackLength, ...
                'MinStepVar',         o.MinStepVar, ...
                'SigmaFixed',         sigma_fixed, ...
                'Verbose',            o.Verbose);
            results.track_type = trackType;

            % Post-fit MinAlpha filter: remove essentially immobile tracks
            if o.MinAlpha > 0
                keep_a = results.alpha >= o.MinAlpha;
                n_removed = sum(~keep_a);
                fields = {'alpha','K','Ka','sigma','track_length','loglik','converged'};
                for fi = 1:numel(fields)
                    if isfield(results, fields{fi})
                        results.(fields{fi}) = results.(fields{fi})(keep_a);
                    end
                end
                results.n_tracks  = sum(keep_a);
                results.n_removed = n_removed;
                if n_removed > 0
                    fprintf('getFBMAlphaDistribution: removed %d tracks with alpha < %.2f\n', ...
                        n_removed, o.MinAlpha);
                end
            else
                results.n_removed = 0;
            end

            if ~isSubset
                obj.FBMAlphaResults    = results;
                obj.IsFBMAlphaComputed = true;
            end
        end

        %% Displacement autocorrelation (DACF) diagnostic
        function [dacf, lags_out] = computeDACF(obj, varargin)
            % COMPUTEDACF  Empirical displacement autocorrelation function.
            %
            %   [dacf, lags] = tc.computeDACF()
            %   [dacf, lags] = tc.computeDACF('MaxLag', 5, 'MinStepVar', 4e-4)
            %
            %   Accumulates the displacement ACF across all tracks.  Near-immobile
            %   tracks (mean squared step < MinStepVar) are excluded to prevent
            %   division-by-near-zero biasing the mean DACF value.
            %
            %   Options:
            %     'TrackType'    - 'culled' (default) | 'raw' | 'relative'
            %     'Condition'    - filter by condition
            %     'MaxLag'       - maximum lag in frames (default: 5)
            %     'MinTrackLength' - minimum track length (default: 10)
            %     'MinStepVar'   - minimum mean squared step (µm²) to include
            %                      track; set to ~4*sigma^2 (default: 0 = no filter)
            %
            %   Output:
            %     dacf     - MaxLag x 1 vector, dacf(k) = DACF at lag k
            %     lags_out - lag indices (1:MaxLag)

            p = inputParser;
            addParameter(p, 'TrackType',    'culled', @(x) ismember(lower(char(x)), {'culled','raw','relative'}));
            addParameter(p, 'Condition',    [],       @(x) ischar(x)||isstring(x)||iscell(x)||isempty(x));
            addParameter(p, 'MaxLag',       5,        @isnumeric);
            addParameter(p, 'MinTrackLength', 10,     @isnumeric);
            addParameter(p, 'MinStepVar',   0,        @isnumeric);
            parse(p, varargin{:});
            o = p.Results;

            trackType = lower(char(o.TrackType));
            [tracks, ~] = obj.getTracksForFBM(trackType, o.Condition);

            maxlag  = round(o.MaxLag);
            Lmin    = round(o.MinTrackLength);
            num     = zeros(maxlag, 1);
            den     = 0;
            n_used  = 0;
            n_skip  = 0;

            for k = 1:numel(tracks)
                tr = tracks{k};
                if size(tr, 1) < Lmin + maxlag
                    continue;
                end
                dr2 = diff(tr(:,1)).^2 + diff(tr(:,2)).^2;  % squared step sizes
                var0 = mean(dr2);

                % Skip near-immobile tracks — their DACF is undefined/noisy
                if o.MinStepVar > 0 && var0 < o.MinStepVar
                    n_skip = n_skip + 1;
                    continue;
                end

                % Use vector components (x and y separately), averaged.
                % Step-magnitude DACF loses sign information and gives near-zero
                % values even for strongly subdiffusive fBM.  Vector-component
                % DACF recovers the correct (2^alpha-2)/2 at lag 1.
                dx = diff(tr(:, 1));
                dy = diff(tr(:, 2));
                v0 = (mean(dx.^2) + mean(dy.^2)) / 2;   % per-component variance
                if v0 < eps
                    n_skip = n_skip + 1;
                    continue;
                end
                for lag = 1:maxlag
                    xc = (mean(dx(1:end-lag) .* dx(1+lag:end)) + ...
                          mean(dy(1:end-lag) .* dy(1+lag:end))) / 2;
                    num(lag) = num(lag) + xc / v0;
                end
                den    = den + 1;
                n_used = n_used + 1;
            end

            if den == 0
                warning('computeDACF: no qualifying tracks found.');
                dacf     = nan(maxlag, 1);
                lags_out = (1:maxlag)';
                return;
            end

            dacf     = num / den;
            lags_out = (1:maxlag)';

            fprintf('computeDACF: %d tracks used, %d skipped (MinStepVar filter)\n', ...
                n_used, n_skip);
            fprintf('  DACF(1) = %.4f  (theoretical for pure fBM: (2^alpha - 2)/2; e.g. -0.34 at alpha=0.4)\n', dacf(1));
        end

        %% Per-RL-state fBM fitting
        function results = getFBMByRLState(obj, varargin)
            % GETFBMBYRLSTATE  Run fBM MLE independently on each RL-classified state.
            %
            %   results = tc.getFBMByRLState()
            %   results = tc.getFBMByRLState('CIMethod','profile', 'MinAlpha',0.1)
            %
            %   Reads tc.RLResults.computed_quantities.classified_tracks and runs
            %   fitFBM_MLE on the tracks belonging to each state.  getRLDecomposition()
            %   must be called first (the RL step already runs a quick fBM fit with
            %   CIMethod='none'; call this method for full profile/MCMC CIs or for
            %   per-track alpha distributions per state).
            %
            %   Options:
            %     'CIMethod'       - 'profile'|'profile+mcmc'|'mcmc'|'hessian'|'none'
            %                        (default: 'profile')
            %     'SubtrackLength' - subtrack length for population MLE (default: 20)
            %     'MinGroupSize'   - minimum tracks per state to attempt fit (default: 50)
            %     'MinStepVar'     - pre-filter within each state (default: 0)
            %     'MinAlpha'       - post-fit alpha floor for per-track results (default: 0)
            %     'PerTrack'       - also run fitFBM_pertracks per state (default: false)
            %     'Verbose'        - print progress (default: false)
            %
            %   Returns a struct array results(s) with one entry per state:
            %     .state        state index
            %     .n_tracks     number of tracks in state
            %     .K, .alpha, .Ka, .sigma   population MLE values
            %     .CI           confidence intervals (from fitFBM_MLE)
            %     .pertracks    per-track struct (if PerTrack=true)
            %     .skipped      true if state was too small to fit

            varargin = unpack_opts(varargin{:});
            p = inputParser;
            addParameter(p, 'CIMethod',          'profile', @ischar);
            addParameter(p, 'SubtrackLength',    20,       @isnumeric);
            addParameter(p, 'MinGroupSize',      50,       @isnumeric);
            addParameter(p, 'MinStepVar',         0,       @isnumeric);
            addParameter(p, 'MinAlpha',           0,       @isnumeric);
            addParameter(p, 'PerTrack',         false,     @islogical);
            addParameter(p, 'MaxSubtrackLength', 20,       @isnumeric);  % matches SubtrackLength by default
            addParameter(p, 'Verbose',          false,     @islogical);
            parse(p, varargin{:});
            o = p.Results;

            if ~obj.IsRLComputed || ~isfield(obj.RLResults, 'computed_quantities')
                error('getFBMByRLState: run tc.getRLDecomposition() first.');
            end

            cq        = obj.RLResults.computed_quantities;
            dt        = obj.RLResults.dt;
            frac      = obj.RLResults.frac;
            rl_type   = obj.RLResults.track_type;
            idx_map   = obj.RLResults.track_index_map;   % filtered→full index translation
            n_states  = numel(cq.classified_tracks);

            % Retrieve the same track population that was used for RL
            switch rl_type
                case 'relative'
                    all_tracks = obj.getAllRelativeTracks();
                case 'raw'
                    all_tracks = obj.getAllRawTracks();
                otherwise
                    all_tracks = obj.getAllCulledTracks();
            end

            results = struct('state', cell(n_states,1), 'n_tracks', cell(n_states,1), ...
                'K', cell(n_states,1), 'alpha', cell(n_states,1), ...
                'Ka', cell(n_states,1), 'sigma', cell(n_states,1), ...
                'CI', cell(n_states,1), 'pertracks', cell(n_states,1), ...
                'skipped', cell(n_states,1));

            for s = 1:n_states
                % classified_tracks{s} indexes into the filtered array passed to RL_HILO.
                % Translate through idx_map to get indices into the full track array.
                rl_idx       = cq.classified_tracks{s};
                idx          = idx_map(rl_idx);
                n_s          = numel(idx);
                results(s).state    = s;
                results(s).n_tracks = n_s;
                results(s).skipped  = false;

                if n_s < o.MinGroupSize
                    fprintf('getFBMByRLState: state %d — %d tracks (< %d), skipped.\n', ...
                        s, n_s, o.MinGroupSize);
                    results(s).skipped = true;
                    continue;
                end

                state_tracks = all_tracks(idx);

                % Optional pre-filter within this state
                if o.MinStepVar > 0
                    mobile = cellfun(@(t) mean(diff(t(:,1)).^2 + diff(t(:,2)).^2) >= o.MinStepVar, state_tracks);
                    state_tracks = state_tracks(mobile);
                    if o.Verbose
                        fprintf('getFBMByRLState: state %d — removed %d stuck tracks\n', s, sum(~mobile));
                    end
                end

                fprintf('getFBMByRLState: state %d (%d tracks) — fitting population fBM MLE...\n', s, numel(state_tracks));
                try
                    [fp, ~, ci] = fitFBM_MLE(state_tracks, dt, ...
                        'ExposureFraction', frac, ...
                        'SubtrackLength',   o.SubtrackLength, ...
                        'CIMethod',         o.CIMethod, ...
                        'Verbose',          o.Verbose);
                    results(s).K     = fp.K;
                    results(s).alpha = fp.alpha;
                    results(s).Ka    = fp.Ka;
                    results(s).sigma = fp.sigma;
                    results(s).CI    = ci;
                    fprintf('  state %d: alpha=%.3f [%.3f, %.3f], K=%.4g, sigma=%.4f µm\n', ...
                        s, fp.alpha, ci.alpha(1), ci.alpha(2), fp.K, fp.sigma);
                catch ME
                    warning('getFBMByRLState: fBM MLE failed for state %d: %s', s, ME.message);
                    results(s).skipped = true;
                end

                % Optional per-track distribution within this state
                if o.PerTrack && ~results(s).skipped
                    pt = fitFBM_pertracks(state_tracks, dt, ...
                        'ExposureFraction',  frac, ...
                        'SigmaFixed',        results(s).sigma, ...
                        'MaxSubtrackLength', o.MaxSubtrackLength, ...
                        'MinStepVar',        o.MinStepVar, ...
                        'Verbose',           o.Verbose);
                    if o.MinAlpha > 0
                        keep_a = pt.alpha >= o.MinAlpha;
                        fields = {'alpha','K','Ka','sigma','track_length','loglik','converged','original_index'};
                        for fi = 1:numel(fields)
                            if isfield(pt, fields{fi})
                                pt.(fields{fi}) = pt.(fields{fi})(keep_a);
                            end
                        end
                        pt.n_tracks  = sum(keep_a);
                        pt.n_removed = sum(~keep_a);
                    end
                    results(s).pertracks = pt;
                    fprintf('  state %d per-track: median alpha=%.3f\n', s, median(pt.alpha, 'omitnan'));
                end
            end
        end

        %% Bootstrap CI (separate from pEM)
        function ciResults = computeBootstrapCI(obj, varargin)
            % COMPUTEBOOTSTRAPCI  Compute confidence intervals on a completed pEM run.
            %
            %   ciResults = tc.computeBootstrapCI()
            %   ciResults = tc.computeBootstrapCI('nBoot', 200, 'parallel', false)
            %
            %   Must be called after getBayesianDiffusivity() has completed for
            %   the full collection (not a condition subset).  Results are stored
            %   in pEMResults.ciresults and included in subsequent calls to
            %   getBayesianDiffusivity().
            %
            %   Name-value options:
            %     'nBoot'          - bootstrap resamples          (default 100)
            %     'covModel'       - 'normal','confined','fBM'    (default 'normal')
            %     'alpha'          - CI level: 1-alpha            (default 0.05)
            %     'nRandomStarts'  - restarts per replicate        (default 5)
            %     'resampleLevel'  - 'parent' or 'subtracks'      (default 'parent')
            %     'parallel'       - use parfor                   (default true)
            %     'verbose'        - print progress               (default true)

            if ~obj.IspEMComputed
                error('TrajectoryCollection:pEMNotComputed', ...
                    'Run getBayesianDiffusivity() on the full collection first.');
            end
            if ~isfield(obj.pEMBootstrapInputs, 'deltaX') || isempty(obj.pEMBootstrapInputs.deltaX)
                error('TrajectoryCollection:noBootstrapInputs', ...
                    ['Bootstrap inputs not found. This can happen if pEM was run on ' ...
                     'a condition subset only. Run getBayesianDiffusivity() on the ' ...
                     'full collection (no Condition argument) first.']);
            end

            p = inputParser;
            addParameter(p, 'nBoot',         100,      @isnumeric);
            addParameter(p, 'covModel',      'normal',  @ischar);
            addParameter(p, 'alpha',          0.05,     @isnumeric);
            addParameter(p, 'nRandomStarts',  5,        @isnumeric);
            addParameter(p, 'resampleLevel', 'parent',  @ischar);
            addParameter(p, 'parallel',       true,     @islogical);
            addParameter(p, 'verbose',        true,     @islogical);
            parse(p, varargin{:});
            o = p.Results;

            deltaX     = obj.pEMBootstrapInputs.deltaX;
            rawResults = obj.pEMBootstrapInputs.rawResults;
            trackInfo  = obj.pEMBootstrapInputs.trackInfo;
            params     = obj.pEMBootstrapInputs.params;

            fprintf('\n=== Bootstrap CI ===\n');
            fprintf('Timestamp  : %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
            fprintf('nBoot      : %d\n', o.nBoot);
            fprintf('covModel   : %s\n', o.covModel);
            fprintf('alpha      : %.3f (%.0f%% CI)\n', o.alpha, 100*(1-o.alpha));
            fprintf('resample   : %s\n', o.resampleLevel);
            fprintf('parallel   : %s\n', ifelse(o.parallel, 'yes', 'no'));

            tic;
            ciResults = pEMv2_bootstrap_CI(deltaX, rawResults, trackInfo, params, ...
                'nBoot',         o.nBoot, ...
                'covModel',      o.covModel, ...
                'alpha',         o.alpha, ...
                'nRandomStarts', o.nRandomStarts, ...
                'resampleLevel', o.resampleLevel, ...
                'parallel',      o.parallel, ...
                'verbose',       o.verbose);
            toc;

            obj.pEMResults.ciresults = ciResults;
            obj.IsBootstrapComputed  = true;
            fprintf('Bootstrap CI stored in pEMResults.ciresults.\n');
            fprintf('===================\n\n');
        end

        %% Lifetime
        function results = getLifetime(obj, varargin)
            if ~obj.IsLifetimeComputed
                p = inputParser;
                addParameter(p, 'rmin',           0.25);
                addParameter(p, 'rmax',           0.4);
                addParameter(p, 'minTrackLength', 4);
                addParameter(p, 'maxTrackLength', 400);
                addParameter(p, 'UseRaw',         false);
                parse(p, varargin{:});

                if p.Results.UseRaw
                    tracks = obj.getAllRawTracks();
                else
                    tracks = obj.getAllCulledTracks();
                end
                if isempty(tracks)
                    error('No tracks available for lifetime computation.');
                end

                tracklen = cellfun('size', tracks, 1);
                X = cell(1, numel(tracks));
                for itrack = 1:numel(tracks)
                    X{itrack} = tracks{itrack}(:, 1:2);   % already in µm
                end

                idx  = find(tracklen > p.Results.minTrackLength & ...
                            tracklen < p.Results.maxTrackLength);
                tmp{1} = X(idx);
                [tracklen, x, f, flo, fup, jump_hist] = ...
                    find_bound_particles_2(tmp, p.Results.rmin, p.Results.rmax);
                prctile(jump_hist, 99)

                obj.LifetimeResults.x           = x;
                obj.LifetimeResults.f           = f;
                obj.LifetimeResults.flo         = flo;
                obj.LifetimeResults.fup         = fup;
                obj.LifetimeResults.TrackLengths = tracklen;
                obj.IsLifetimeComputed = true;
            end
            results = obj.LifetimeResults;
        end

        function obj = recomputeLifetime(obj, varargin)
            obj.IsLifetimeComputed = false;
            obj.LifetimeResults    = struct();
            obj.getLifetime(varargin{:});
            fprintf('Lifetime recomputed.\n');
        end

        function obj = recomputeRL(obj, varargin)
            obj.IsRLComputed = false;
            obj.RLResults    = struct();
            obj.getRLDecomposition(varargin{:});
            fprintf('RL decomposition recomputed.\n');
        end

        function results = recomputeMSD(obj, varargin)
            results = obj.getMSD(varargin{:}, 'ForceRecompute', true);
            fprintf('MSD recomputed.\n');
        end

        %% Visualization
        function fig = visualizeSubsets(obj, varargin)
            p = inputParser;
            addParameter(p, 'Condition', '');
            addParameter(p, 'FileID',    '');
            addParameter(p, 'NumTracks', inf);
            addParameter(p, 'UseRaw',    false);
            parse(p, varargin{:});

            idx = true(size(obj.Wrappers));
            if ~isempty(p.Results.Condition)
                idx = idx & strcmp(obj.Metadata.Condition, p.Results.Condition);
            end
            if ~isempty(p.Results.FileID)
                idx = idx & strcmp(obj.Metadata.FileID, p.Results.FileID);
            end
            subWrappers = obj.Wrappers(idx);

            fig    = figure;
            hold on;
            colors = lines(length(subWrappers));
            for w = 1:length(subWrappers)
                tracks = ifelse(p.Results.UseRaw, subWrappers{w}.getRawTracks(), ...
                    subWrappers{w}.getCulledTracks());
                if length(tracks) > p.Results.NumTracks
                    tracks = tracks(randperm(length(tracks), p.Results.NumTracks));
                end
                for t = 1:length(tracks)
                    plot(tracks{t}(:,1), tracks{t}(:,2), 'Color', colors(w,:));
                end
                if ~isempty(subWrappers{w}.getNucleusMask())
                    contour(subWrappers{w}.getNucleusMask(), [0.5 0.5], 'k--');
                end
            end
            title('Trajectory Subsets');
            xlabel('X'); ylabel('Y');
            hold off;
        end

        function figs = plotMSD(obj, varargin)
            % PLOTMSD  Plot ensemble MSD with bootstrap CI and median fBM fit.
            %
            %   figs = tc.plotMSD()
            %   figs = tc.plotMSD('Color',  [0 0.4 0.8])
            %   figs = tc.plotMSD('Title',  'H2B Control')
            %   figs = tc.plotMSD('ShowBootstraps', true)   % overlay individual fits
            %
            % Returns a figure handle.  Call tc.getMSD() first.

            if ~obj.IsMSDComputed
                error('plotMSD: run tc.getMSD() first.');
            end

            varargin = unpack_opts(varargin{:});
            p = inputParser;
            addParameter(p, 'Color',          [0 0.4 0.8], @isnumeric);
            addParameter(p, 'Title',          '',          @(x) ischar(x)||isstring(x));
            addParameter(p, 'ShowBootstraps', false,       @islogical);
            parse(p, varargin{:});

            r      = obj.MSDResults;
            tau    = r.lag_axis;
            emsd   = r.ensemble_mean;
            bmeans = r.boot_means;
            bfits  = r.boot_fits;
            dt     = r.dt;
            frac   = r.frac;
            lags   = r.lag_frames;
            clr    = p.Results.Color;

            ci_lo = prctile(bmeans, 2.5,  2);
            ci_hi = prctile(bmeans, 97.5, 2);

            figs = figure('Name', 'Ensemble MSD');
            patch([tau; flipud(tau)], [ci_lo; flipud(ci_hi)], clr, ...
                  'EdgeColor', 'none', 'FaceAlpha', 0.25);
            hold on;

            % Individual bootstrap fit curves (optional)
            if p.Results.ShowBootstraps
                for b = 1:size(bfits, 1)
                    fp = bfits(b, :);
                    if any(isnan(fp)), continue; end
                    fc = msd_model_curve(fp, lags, dt, frac);
                    plot(tau, fc, 'Color', [clr 0.08], 'LineWidth', 0.5);
                end
            end

            % Ensemble mean
            plot(tau, emsd, '-', 'Color', clr, 'LineWidth', 2.5);

            % Median fit curve
            valid = bfits(~any(isnan(bfits), 2), :);
            if ~isempty(valid)
                med_fp   = median(valid, 1);
                fit_crv  = msd_model_curve(med_fp, lags, dt, frac);
                plot(tau, fit_crv, '--', 'Color', clr * 0.6, 'LineWidth', 1.5);
                med_alpha = med_fp(3);
                ci_alpha  = prctile(valid(:,3), [2.5 97.5]);
                legend({'95% CI', 'Ensemble mean', ...
                        sprintf('fBM fit  \\alpha=%.2f [%.2f–%.2f]', ...
                            med_alpha, ci_alpha(1), ci_alpha(2))}, ...
                       'Location', 'northwest', 'Box', 'off');
            else
                legend({'95% CI', 'Ensemble mean'}, 'Location', 'northwest', 'Box', 'off');
            end

            set(gca, 'XScale', 'log', 'YScale', 'log', ...
                     'FontSize', 18, 'FontWeight', 'bold', 'LineWidth', 1.5);
            xlabel('\tau (s)');
            ylabel('MSD (\mum^2)');
            ttl = char(p.Results.Title);
            if isempty(ttl)
                ttl = sprintf('Ensemble MSD  (n=%d tracks)', r.n_tracks);
            end
            title(ttl);
            axis square; box off;
        end

        function figs = plotFBMFit(obj, varargin)
            % PLOTFBMFIT  Overlay the Toeplitz MLE fBM curve on the empirical
            %             ensemble MSD.
            %
            %   figs = tc.plotFBMFit('FBMResults', fbm)   % pass result from getFBMParameters
            %   figs = tc.plotFBMFit('Color', [0 0.4 0.8])
            %   figs = tc.plotFBMFit('Title', 'H2B Control')
            %   figs = tc.plotFBMFit('ShowBootstrap', true)
            %
            % Requires tc.getMSD() to have been run.  Pass 'FBMResults' explicitly
            % when getFBMParameters was called with TrackType='relative' or a
            % Condition filter (those results are not cached on the object):
            %
            %   fbm = tc.getFBMParameters('TrackType','relative', ...);
            %   tc.plotFBMFit('FBMResults', fbm)
            %
            % The Toeplitz MLE curve uses the exact Backlund exposure-time kernel:
            %   MSD_2d(n) = 4K * (psi(n*dt, Te, alpha) - psi(0, Te, alpha)) + 4*sigma^2

            varargin = unpack_opts(varargin{:});
            if ~obj.IsMSDComputed
                error('plotFBMFit: run tc.getMSD() first.');
            end

            p = inputParser;
            addParameter(p, 'FBMResults',    [],          @(x) isstruct(x)||isempty(x));
            addParameter(p, 'Color',         [0 0.4 0.8], @isnumeric);
            addParameter(p, 'Title',         '',          @(x) ischar(x)||isstring(x));
            addParameter(p, 'ShowBootstrap', true,        @islogical);
            parse(p, varargin{:});

            % ── Resolve FBM results ───────────────────────────────────────────
            fb = p.Results.FBMResults;
            if isempty(fb)
                if ~obj.IsFBMComputed
                    error(['plotFBMFit: no FBM results available. Either run ' ...
                           'tc.getFBMParameters() (for culled tracks) or pass ' ...
                           '''FBMResults'' explicitly when using TrackType=''relative''.']);
                end
                fb = obj.FBMResults;
            end

            % ── MSD data ──────────────────────────────────────────────────────
            r      = obj.MSDResults;
            tau    = r.lag_axis;       % lag times in seconds
            emsd   = r.ensemble_mean;
            bmeans = r.boot_means;
            bfits  = r.boot_fits;
            dt     = r.dt;
            frac   = r.frac;
            lags   = r.lag_frames;     % integer frame lags (1, 2, …, MaxLag)
            clr    = p.Results.Color;

            ci_lo = prctile(bmeans, 2.5,  2);
            ci_hi = prctile(bmeans, 97.5, 2);

            % ── Toeplitz MLE parameters ───────────────────────────────────────
            K     = fb.K;
            alpha = fb.alpha;
            sigma = fb.sigma;
            % Use the ExposureFraction from the FBM fit, NOT the MSD frac.
            if isfield(fb, 'frac')
                Te = fb.frac * dt;
            else
                Te = frac * dt;
            end

            % Backlund psi function (exact exposure-time kernel)
            psi = @(tau_s) ...
                ((tau_s + Te).^(alpha+2) + abs(tau_s - Te).^(alpha+2) ...
                 - 2*tau_s.^(alpha+2)) ./ (Te^2 * (alpha+1) * (alpha+2));

            psi_0   = 2 * Te^alpha / ((alpha+1) * (alpha+2));
            tau_s   = lags * dt;       % lag times in seconds (column vector)
            msd_mle = 4*K * (psi(tau_s) - psi_0) + 4*sigma^2;

            % ── Plot ──────────────────────────────────────────────────────────
            figs = figure('Name', 'Ensemble MSD with Toeplitz MLE fit');

            % Bootstrap CI shading (from MSD bootstrap)
            patch([tau; flipud(tau)], [ci_lo; flipud(ci_hi)], clr, ...
                  'EdgeColor', 'none', 'FaceAlpha', 0.25);
            hold on;

            % MSD least-squares fit curve (median over bootstrap, dashed)
            if p.Results.ShowBootstrap
                valid = bfits(~any(isnan(bfits), 2), :);
                if ~isempty(valid)
                    med_fp  = median(valid, 1);
                    fit_crv = msd_model_curve(med_fp, lags, dt, frac);
                    plot(tau, fit_crv, '--', 'Color', clr * 0.6, 'LineWidth', 1.5, ...
                         'DisplayName', sprintf('MSD fit \\alpha=%.2f', med_fp(3)));
                end
            end

            % Empirical ensemble mean
            plot(tau, emsd, '-', 'Color', clr, 'LineWidth', 2.5, ...
                 'DisplayName', 'Ensemble mean');

            % Toeplitz MLE curve (black dashed)
            plot(tau, msd_mle, 'k--', 'LineWidth', 2, ...
                 'DisplayName', sprintf('Toeplitz MLE \\alpha=%.3f, K=%.4f, \\sigma=%.4f', ...
                                        alpha, K, sigma));

            set(gca, 'XScale', 'log', 'YScale', 'log', ...
                     'FontSize', 18, 'FontWeight', 'bold', 'LineWidth', 1.5);
            xlabel('\tau (s)');
            ylabel('MSD (\mum^2)');
            ttl = char(p.Results.Title);
            if isempty(ttl)
                ttl = sprintf('MSD + Toeplitz MLE  (n=%d tracks)', r.n_tracks);
            end
            title(ttl);
            legend('Location', 'northwest', 'Box', 'off');
            axis square; box off;
        end

        function figs = plotRLDecomposition(obj, varargin)
            % PLOTRLDECOMPOSITION  Reproduce the standard RL analysis figures
            %                      from stored results without re-running RL_HILO.
            %
            %   figs = tc.plotRLDecomposition()
            %   figs = tc.plotRLDecomposition('Title', 'H2B Control')
            %
            % Produces three figures:
            %   Fig 1 – van Hove correlation G(r,τ): data + RL fit
            %   Fig 2 – P(MSD): Richardson-Lucy deconvolved MSD distribution
            %   Fig 3 – Per-state ensemble MSD with 99% CI error bars
            %
            % Call tc.getRLDecomposition() first.

            if ~obj.IsRLComputed
                error('plotRLDecomposition: run tc.getRLDecomposition() first.');
            end

            varargin = unpack_opts(varargin{:});
            p = inputParser;
            addParameter(p, 'Title', '', @(x) ischar(x)||isstring(x));
            parse(p, varargin{:});

            cq      = obj.RLResults.computed_quantities;
            ttl     = char(p.Results.Title);
            if isempty(ttl), ttl = cq.ident; end

            dt      = obj.getFrameIntervalForCondition([]);
            n_lags  = size(cq.msd, 1);
            tau_rl  = (1:n_lags)' * dt;
            n_states = numel(cq.classified_tracks);
            colors   = lines(n_states);

            % --- Fig 1: van Hove correlation ---------------------------------
            figs(1) = figure('Name', 'van Hove correlation');
            plot(cq.x, cq.vanHove, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'none');
            hold on;
            plot(cq.x, cq.Gs, 'r-', 'LineWidth', 2);
            set(gca, 'YScale', 'log', 'FontSize', 24, 'FontWeight', 'bold', 'LineWidth', 1);
            xlabel('r  (\mum)');
            ylabel('G(r,\tau)');
            axis([0 1 1e-4 100]);
            axis square; box off;
            title(ttl);
            legend({'Data', 'RL fit'}, 'Location', 'northeast', 'Box', 'off');

            % --- Fig 2: P(MSD) -----------------------------------------------
            figs(2) = figure('Name', 'P(MSD)');
            semilogx(obj.RLResults.M, cq.P1norm, 'k-', 'LineWidth', 2);
            set(gca, 'FontSize', 24, 'FontWeight', 'bold', 'LineWidth', 1);
            xlabel('MSD  (\mum^2)');
            ylabel('P(MSD)');
            axis([1e-3 2e-1 0 max(cq.P1norm) * 1.1]);
            title(ttl);
            box off;

            % --- Fig 3: per-state ensemble MSD -------------------------------
            figs(3) = figure('Name', 'Per-state MSD');
            hold on;
            for s = 1:n_states
                n_s = numel(cq.classified_tracks{s});
                if n_s == 0, continue; end
                errorbar(tau_rl, cq.msd(:, s), cq.msderr(:, s), ...
                         '-o', 'Color', colors(s, :), 'LineWidth', 1.5, ...
                         'MarkerSize', 4, 'CapSize', 3, ...
                         'DisplayName', sprintf('State %d  (n=%d)', s, n_s));
            end
            set(gca, 'XScale', 'log', 'YScale', 'log', ...
                     'FontSize', 18, 'FontWeight', 'bold', 'LineWidth', 1.5);
            xlabel('\tau (s)');
            ylabel('MSD  (\mum^2)');
            title([ttl '  —  per-state MSD']);
            legend('show', 'Location', 'northwest', 'Box', 'off');
            axis square; box off;
        end

        function globalfig = plotpEMstats(obj)
            if ~obj.IspEMComputed
                warning('First compute pEM');
                return;
            end
            proportions  = obj.pEMResults.proportions;
            pp           = obj.pEMResults.pEMTable.posteriorProb{1};
            max_pp       = max(pp, [], 2);
            sort_pp      = sort(pp, 2, 'descend');
            delta_pp     = sort_pp(:,1) - sort_pp(:,2);
            idx          = find(max_pp > 0.5 & delta_pp > 0.2);
            optimalState = obj.pEMResults.pEMTable.optimalState{1};
            optimalState(setdiff(1:length(optimalState), idx)) = 9;
            actual_num   = length(optimalState) - length(find(optimalState == 9));
            kept_states  = true(height(proportions), 1);
            for ii = 1:height(proportions)
                if length(find(optimalState==ii)) / actual_num < 0.025
                    optimalState(find(optimalState==ii)) = 9;
                    kept_states(ii) = false;
                else
                    sprintf('state %d proportion %f', ii, length(find(optimalState==ii)) / actual_num);
                end
            end
            splitLength       = 7;
            colors            = [[0.7 0.3 0.3];[1 0 0];[0 0 1];[0 0.8 0.2];[0 1 0];[0.2 0.2 0.2]];
            msd_all_states    = zeros(splitLength-1, height(proportions));
            var_all_states    = zeros(splitLength-1, height(proportions));
            for ii = 1:length(kept_states)
                if kept_states(ii)
                    idx2 = find(optimalState == ii);
                    dum  = obj.pEMResults.pEMTable.splitX{1}(idx2);
                    ensemble_msd  = zeros(splitLength-1, length(idx2));
                    ensemble_lags = zeros(splitLength-1, length(idx2));
                    for jj = 1:length(dum)
                        [msd, num_lags_tracks, ~] = simple_msd(dum{jj}(:,1:2), (1:7)');
                        ensemble_msd(:,jj)  = msd;
                        ensemble_lags(:,jj) = num_lags_tracks;
                    end
                    msd_all_states(:,ii) = mean(ensemble_msd ./ ensemble_lags, 2);
                    var_all_states(:,ii) = std(ensemble_msd ./ ensemble_lags, 0, 2) ./ ...
                                           sqrt(sum(ensemble_lags, 2));
                end
            end
            globalfig(1) = figure('Position', [1 1 0.5 1] .* get(0,'Screensize'), 'Visible','on');
            dt = obj.getFrameIntervalForCondition([]);
            for ii = 1:length(kept_states)
                if kept_states(ii)
                    errorbar(dt*(0:6), [0; msd_all_states(:,ii)], [0; var_all_states(:,ii)], ...
                             's', 'Color', colors(ii,:), 'markersize', 8);
                    hold on;
                    plot(dt*(0:6), [0; msd_all_states(:,ii)], 'Color', colors(ii,:), 'linewidth', 1);
                end
            end
            axis square; box off;
            set(gca, 'linewidth', 2, 'fontsize', 24, 'fontweight', 'bold');
        end

        function globalfig = plotLifetime(obj)
            x   = obj.LifetimeResults.x;
            f   = obj.LifetimeResults.f;
            flo = obj.LifetimeResults.flo;
            fup = obj.LifetimeResults.fup;
        end

        %% Metadata / query helpers
        function obj = setCondition(obj, fileID, newCondition)
            idx = strcmp(obj.Metadata.FileID, fileID);
            if any(idx)
                obj.Metadata.Condition(idx) = newCondition;
            else
                warning('FileID not found.');
            end
        end

        function tracklen = getTrackLengthDistribution(obj)
            tracks   = obj.getAllCulledTracks();
            tracklen = cellfun('size', tracks, 1);
        end

        function tracks = getTracksByCondition(obj, condition, varargin)
            % GETTRACKSBYCONDITION  Return culled (or raw) tracks for a condition.
            %
            %   tracks = tc.getTracksByCondition('Control')
            %   tracks = tc.getTracksByCondition('Treated', 'UseRaw', true)

            p = inputParser;
            addRequired(p, 'condition');
            addParameter(p, 'UseRaw', false, @islogical);
            parse(p, condition, varargin{:});

            if ischar(condition) || isstring(condition)
                condition = {char(condition)};
            end

            matchIdx = false(size(obj.Metadata.Condition));
            for c = condition(:)'
                matchIdx = matchIdx | strcmp(obj.Metadata.Condition, c);
            end

            if ~any(matchIdx)
                warning('No files found for condition(s): %s', strjoin(condition, ', '));
                tracks = {};
                return;
            end

            matchingWrappers = obj.Wrappers(matchIdx);
            tracks = cell(0,1);

            for i = 1:length(matchingWrappers)
                tw = matchingWrappers{i};
                if p.Results.UseRaw
                    thisTracks = tw.getRawTracks();
                else
                    thisTracks = tw.getCulledTracks();
                end
                if ~isempty(thisTracks)
                    if ~iscolumn(thisTracks), thisTracks = thisTracks(:); end
                    tracks = [tracks; thisTracks]; %#ok<AGROW>
                end
            end

            fprintf('Returned %d tracks for condition(s): %s\n', ...
                length(tracks), strjoin(condition, ', '));
        end

        function tracks = getRelativeTracksByCondition(obj, condition)
            % GETRELATIVETRACKSBYCONDITION  Return relative tracks for a condition.
            if ischar(condition) || isstring(condition)
                condition = {char(condition)};
            end
            matchIdx = false(size(obj.Metadata.Condition));
            for c = condition(:)'
                matchIdx = matchIdx | strcmp(obj.Metadata.Condition, c);
            end
            if ~any(matchIdx)
                warning('No files found for condition(s): %s', strjoin(condition, ', '));
                tracks = {};
                return;
            end
            tracks = cell(0, 1);
            for tw = obj.Wrappers(matchIdx)
                rt = tw{1}.getRelativeTracks();
                if ~isempty(rt)
                    if ~iscolumn(rt), rt = rt(:); end
                    tracks = [tracks; rt]; %#ok<AGROW>
                end
            end
        end

        function summary(obj)
            fprintf('\n=== TrajectoryCollection Summary ===\n');
            fprintf('Number of files/wrappers: %d\n', length(obj.Wrappers));

            if ~isempty(obj.Wrappers)
                rawCounts    = cellfun(@(tw) tw.getNumRawTracks(),    obj.Wrappers);
                culledCounts = cellfun(@(tw) tw.getNumCulledTracks(), obj.Wrappers);
                relCounts    = cellfun(@(tw) numel(tw.getRelativeTracks()), obj.Wrappers);
                srcTypes     = cellfun(@(tw) tw.getSourceType(),       obj.Wrappers, 'UniformOutput', false);
                fis          = cellfun(@(tw) tw.getFrameInterval(),    obj.Wrappers);

                fprintf('Total raw tracks     : %d\n', sum(rawCounts));
                fprintf('Total culled tracks  : %d\n', sum(culledCounts));
                fprintf('Total relative tracks: %d\n', sum(relCounts));

                fprintf('\nPer-wrapper details:\n');
                for w = 1:length(obj.Wrappers)
                    fiStr = ifelse(isnan(fis(w)), 'NOT SET', sprintf('%.4f s', fis(w)));
                    fprintf('  [%s] source=%-4s  culled=%4d  px=%.4f  dt=%s\n', ...
                        char(obj.Metadata.FileID(w)), srcTypes{w}, ...
                        culledCounts(w), obj.Wrappers{w}.getPixelSize(), fiStr);
                end
            else
                fprintf('Total raw tracks     : 0\n');
                fprintf('Total culled tracks  : 0\n');
                fprintf('Total relative tracks: 0\n');
            end

            uniqueConditions = unique(obj.Metadata.Condition);
            if ~isempty(uniqueConditions) && ~all(strcmp(uniqueConditions, ''))
                fprintf('Conditions           : %s\n', strjoin(uniqueConditions, ', '));
            else
                fprintf('Conditions           : (none specified)\n');
            end

            fprintf('RL decomp computed   : %s\n', ifelse(obj.IsRLComputed,        'Yes', 'No'));
            fprintf('MSD computed         : %s\n', ifelse(obj.IsMSDComputed,       'Yes', 'No'));
            fprintf('Bayesian computed    : %s\n', ifelse(obj.IspEMComputed,       'Yes', 'No'));
            fprintf('Bootstrap CI computed: %s\n', ifelse(obj.IsBootstrapComputed, 'Yes', 'No'));
            fprintf('====================================\n\n');
        end

        function exportToCSV(obj, outputFile, varargin)
            p = inputParser;
            addParameter(p, 'UseRaw', false);
            parse(p, varargin{:});

            tracks  = ifelse(p.Results.UseRaw, obj.getAllRawTracks(), obj.getAllCulledTracks());
            allData = [];
            for i = 1:length(tracks)
                meta    = repmat(obj.Metadata(i,:), size(tracks{i},1), 1);
                allData = [allData; [tracks{i}, meta.FileID, meta.Condition]]; %#ok<AGROW>
            end
            writetable(array2table(allData), outputFile);
        end

        function tw = getWrapperForTrack(obj, trackIdx, varargin)
            p = inputParser;
            addParameter(p, 'UseRaw', false, @islogical);
            parse(p, varargin{:});

            tracks     = ifelse(p.Results.UseRaw, obj.getAllRawTracks(), obj.getAllCulledTracks());
            cumulative = cumsum(cellfun(@(w) w.getNumCulledTracks(), obj.Wrappers));
            if p.Results.UseRaw
                cumulative = cumsum(cellfun(@(w) w.getNumRawTracks(), obj.Wrappers));
            end

            wrapperIdx = find(cumulative >= trackIdx, 1, 'first');
            if isempty(wrapperIdx)
                error('Track index %d out of range (max %d).', trackIdx, length(tracks));
            end
            tw = obj.Wrappers{wrapperIdx};
        end

        function singleTrack = getTrack(obj, trackIdx, varargin)
            p = inputParser;
            addParameter(p, 'UseRaw', false);
            parse(p, varargin{:});

            allTracks = ifelse(p.Results.UseRaw, obj.getAllRawTracks(), obj.getAllCulledTracks());
            if trackIdx < 1 || trackIdx > length(allTracks)
                error('Track index %d out of range (1 to %d).', trackIdx, length(allTracks));
            end
            singleTrack = allTracks{trackIdx};
        end

        %% Save
        function save(obj, filename)
            % SAVE  Persist the TrajectoryCollection to a .mat file.
            %
            %   tc.save('results.mat')
            %
            %   Reload with:
            %     data = load('results.mat');  tc = data.tc;
            %     tc.reloadTracks()            % optional: rebuild track data from CSVs
            %
            %   Track matrices (RawTracks / RelativeTracks) are automatically
            %   stripped from each Wrapper at serialization time via saveobj —
            %   the live in-memory object is NOT modified.  All analysis results
            %   (pEM, bootstrap CI, MSD, lifetime) and enough metadata to
            %   reconstruct tracks (file paths, ROIs, Parameters, PixelSize,
            %   FrameInterval) are preserved.
            %
            %   Typical file sizes:
            %     Without tracks : ~100-400 MB  (dominated by pEMBootstrapInputs)
            %     tc.reloadTracks() : 5-10 min to rebuild from original CSVs

            if ~endsWith(filename, '.mat')
                filename = [filename '.mat'];
            end

            % Clear TC-level track caches — they are derived from Wrapper tracks
            % and will be stale after a save/load cycle anyway.
            obj.AllRawTracks           = {};
            obj.AllCulledTracks        = {};
            obj.AllRelativeTracks      = {};
            obj.IsAllRawCollected      = false;
            obj.IsAllCulledCollected   = false;
            obj.IsAllRelativeCollected = false;

            tc = obj; %#ok<NASGU>
            saveTCToFile(filename, tc);
            fprintf('Saved TrajectoryCollection (%d wrappers) to %s\n', ...
                numel(obj.Wrappers), filename);
        end

        %% Reload tracks from original CSV files
        function reloadTracks(obj, varargin)
            % RELOADTRACKS  Rebuild RawTracks in each Wrapper from original CSV files.
            %
            %   tc.reloadTracks()
            %   tc.reloadTracks('OldRoot', '/old/mount/', 'NewRoot', '/new/mount/')
            %
            %   Call this after loading a saved TrajectoryCollection to restore
            %   track data.  Analysis results (pEM, MSD, etc.) are unaffected.
            %
            %   If CSV files were added with absolute paths (recommended), no
            %   arguments are needed — reloadTracks() finds them automatically.
            %
            %   Name-value options:
            %     'OldRoot' / 'NewRoot' - pair that replaces a path prefix when
            %                    data has moved (e.g. different mount point).
            %                    Both must be supplied together.
            %                    E.g.: 'OldRoot','/lab/server/', 'NewRoot','/Volumes/Lab/'
            %                    preserves subdirectory structure under the new root.

            p = inputParser;
            addParameter(p, 'OldRoot', '', @ischar);
            addParameter(p, 'NewRoot', '', @ischar);
            parse(p, varargin{:});
            oldRoot = p.Results.OldRoot;
            newRoot = p.Results.NewRoot;

            nReloaded = 0;
            for i = 1:numel(obj.Wrappers)
                tw = obj.Wrappers{i};

                if ~isempty(tw.getRawTracks())
                    continue;   % already has tracks
                end

                if ~strcmp(tw.SourceType, 'csv')
                    warning('TrajectoryCollection:reloadTracks', ...
                        'Wrapper %d is SMD-sourced; tracks cannot be reloaded from CSV.', i);
                    continue;
                end

                fname = tw.FileName;
                if ~isempty(oldRoot) && ~isempty(newRoot) && startsWith(fname, oldRoot)
                    fname = fullfile(newRoot, fname(numel(oldRoot)+1:end));
                end

                if isempty(fname) || ~isfile(fname)
                    warning('TrajectoryCollection:reloadTracks', ...
                        'File not found for wrapper %d: %s', i, fname);
                    continue;
                end

                wasRelative = tw.IsRelative;
                tw.readData(fname);
                % ROIs were preserved by saveobj/loadobj — no need to re-run loadNucleusMask
                if tw.IsCulled || ~isempty(tw.ROIs)
                    tw.cull_tracks();
                    if wasRelative
                        tw.remove_relative_motion();
                    end
                end
                nReloaded = nReloaded + 1;
            end

            obj.invalidateCollections();
            fprintf('Reloaded tracks for %d/%d wrappers.\n', nReloaded, numel(obj.Wrappers));
        end

    end % public methods

    methods (Access = private)

        function collectAllTracks(obj)
            nWrappers = length(obj.Wrappers);
            if nWrappers == 0
                obj.AllRawTracks       = {};
                obj.AllCulledTracks    = {};
                obj.AllRelativeTracks  = {};
                obj.IsAllRawCollected      = true;
                obj.IsAllCulledCollected   = true;
                obj.IsAllRelativeCollected = true;
                return;
            end

            raw      = cell(0,1);
            culled   = cell(0,1);
            relative = cell(0,1);

            for i = 1:nWrappers
                tw           = obj.Wrappers{i};
                thisRaw      = tw.getRawTracks();
                thisCulled   = tw.getCulledTracks();
                thisRelative = tw.getRelativeTracks();

                if ~isempty(thisRaw)      && ~iscolumn(thisRaw),      thisRaw      = thisRaw(:);      end
                if ~isempty(thisCulled)   && ~iscolumn(thisCulled),   thisCulled   = thisCulled(:);   end
                if ~isempty(thisRelative) && ~iscolumn(thisRelative),  thisRelative = thisRelative(:); end

                raw      = [raw;      thisRaw];      %#ok<AGROW>
                culled   = [culled;   thisCulled];   %#ok<AGROW>
                relative = [relative; thisRelative]; %#ok<AGROW>
            end

            obj.AllRawTracks       = raw;
            obj.AllCulledTracks    = culled;
            obj.AllRelativeTracks  = relative;
            obj.IsAllRawCollected      = true;
            obj.IsAllCulledCollected   = true;
            obj.IsAllRelativeCollected = true;

            fprintf('Collected %d raw, %d culled, and %d relative tracks across %d FOVs.\n', ...
                length(raw), length(culled), length(relative), nWrappers);
        end

        function dt = getFrameIntervalForCondition(obj, condition)
            % Returns the frame interval (seconds) for the requested condition.
            % Errors if any matching wrapper has FrameInterval = NaN.
            % Warns if wrappers disagree on the value.

            if isempty(condition)
                wrappers = obj.Wrappers;
            else
                if ischar(condition) || isstring(condition)
                    condition = {char(condition)};
                end
                matchIdx = false(size(obj.Metadata.Condition));
                for c = condition(:)'
                    matchIdx = matchIdx | strcmp(obj.Metadata.Condition, c);
                end
                wrappers = obj.Wrappers(matchIdx);
            end

            if isempty(wrappers)
                error('TrajectoryCollection:noWrappers', ...
                      'No wrappers found for the specified condition.');
            end

            fis = cellfun(@(tw) tw.getFrameInterval(), wrappers);

            if any(isnan(fis))
                error('TrajectoryCollection:frameIntervalNotSet', ...
                    ['FrameInterval is NaN for %d wrapper(s). ' ...
                     'Call tw.readParams(tomlFile) or set FrameInterval explicitly ' ...
                     'before running analysis.'], sum(isnan(fis)));
            end

            uniqueFIs = unique(fis);
            if numel(uniqueFIs) > 1
                warning('TrajectoryCollection:inconsistentFrameInterval', ...
                    'Wrappers have %d different FrameInterval values: [%s]. Using median = %.4f s.', ...
                    numel(uniqueFIs), num2str(uniqueFIs(:)', '%.4f '), median(fis));
            end

            dt = median(fis);
        end

        function invalidateCollections(obj)
            obj.IsAllRawCollected      = false;
            obj.IsAllCulledCollected   = false;
            obj.IsAllRelativeCollected = false;
            obj.IsRLComputed           = false;
            obj.IsMSDComputed          = false;
            obj.IspEMComputed          = false;
            obj.IsBootstrapComputed    = false;
            obj.AllRawTracks       = {};
            obj.AllCulledTracks    = {};
            obj.AllRelativeTracks  = {};
            obj.RLResults              = struct();
            obj.MSDResults             = struct();
            obj.pEMResults             = struct();
            obj.pEMBootstrapInputs     = struct();
            obj.FBMResults             = struct();
            obj.FBMAlphaResults        = struct();
            obj.IsFBMComputed          = false;
            obj.IsFBMAlphaComputed     = false;
        end

        function [tracks, dt] = getTracksForFBM(obj, trackType, condition)
            % GETTRACKSFORFBM  Shared track-getter for fBM methods.
            %   trackType : 'culled' | 'raw' | 'relative'
            %   condition : [] for all, or condition string/cell
            dt = obj.getFrameIntervalForCondition(condition);
            hasCondition = ~isempty(condition);
            switch trackType
                case 'culled'
                    if hasCondition
                        tracks = obj.getTracksByCondition(condition);
                    else
                        tracks = obj.getAllCulledTracks();
                    end
                case 'raw'
                    if hasCondition
                        tracks = obj.getTracksByCondition(condition, 'UseRaw', true);
                    else
                        tracks = obj.getAllRawTracks();
                    end
                case 'relative'
                    if hasCondition
                        tracks = obj.getRelativeTracksByCondition(condition);
                    else
                        tracks = obj.getAllRelativeTracks();
                    end
                otherwise
                    error('getTracksForFBM: unknown TrackType ''%s''.', trackType);
            end
        end

        function addGlobalParameters(obj, params)
            for w = obj.Wrappers
                w{1}.Parameters = mergeStructs(w{1}.Parameters, params);
            end
        end

    end % private methods
end

% -------------------------------------------------------------------------
% Module-level helpers
% -------------------------------------------------------------------------
function str = ifelse(cond, t, f)
    if cond, str = t; else, str = f; end
end

function varargin = unpack_opts(varargin)
% UNPACK_OPTS  Allow methods to accept a struct in place of name-value pairs.
%   If the first (and only) argument is a struct, expand it to name-value pairs.
%   Usage: varargin = unpack_opts(varargin{:});
if numel(varargin) == 1 && isstruct(varargin{1})
    s = varargin{1};
    f = fieldnames(s);
    varargin = reshape([f'; struct2cell(s)'], 1, []);
end
end

function tracks = filter_mobile(tracks, min_mean_sq_step)
% FILTER_MOBILE  Remove tracks whose mean one-step squared displacement is
%   below min_mean_sq_step (µm²).  Typical threshold: 4*sigma_fixed^2 where
%   sigma_fixed is the localization precision measured from immobile particles.
if min_mean_sq_step <= 0
    return;
end
mobile = cellfun(@(t) mean(diff(t(:,1)).^2 + diff(t(:,2)).^2) >= min_mean_sq_step, tracks);
tracks = tracks(mobile);
end

function curve = msd_model_curve(fp, lags, dt, frac)
% Reconstruct fBM MSD model values from fit parameters.
%   fp = [G, sigma_sq, alpha]; lags = integer lag indices.
G = fp(1); sig2 = fp(2); alpha = fp(3);
b     = (abs(1 + frac./lags).^(2+alpha) + abs(1-frac./lags).^(2+alpha) - 2) ./ (frac./lags).^2;
curve = G / ((1+alpha)*(2+alpha)) .* ((dt.*lags).^alpha .* b - 2*(frac*dt)^alpha) + 2*sig2;
end

function s = mergeStructs(s1, s2)
    s      = s1;
    fields = fieldnames(s2);
    for i = 1:length(fields)
        s.(fields{i}) = s2.(fields{i});
    end
end
