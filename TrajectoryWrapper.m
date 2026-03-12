classdef TrajectoryWrapper < handle
    % TRAJECTORYWRAPPER - Modular loader and filter for single-molecule trajectories
    %
    % Supports two input modes:
    %   CSV path  : tw = TrajectoryWrapper(); tw.readData('locs.csv');
    %   SMD path  : tw = TrajectoryWrapper.fromSMD(smd);
    %
    % Common workflow (CSV):
    %   tw = TrajectoryWrapper('PixelSize', 0.165, 'FrameInterval', 0.02);
    %   tw.readData('locs.csv');
    %   tw.loadNucleusMask();
    %   tw.cull_tracks();
    %   tw.remove_relative_motion();
    %   tracks = tw.getCulledTracks();
    %
    % Common workflow (SMD):
    %   smd.localize(); smd.track(); smd.get_roi(); smd.cull_tracks();
    %   tw = TrajectoryWrapper.fromSMD(smd);

    properties (SetAccess = private)
        RawTracks        % Cell array: all trajectories (Nx4: x,y,frame,SNR)
        CulledTracks     % Cell array: filtered trajectories (Nx5: x,y,frame,SNR,ROI_ID)
        RelativeTracks   % Cell array: tracks with global cell motion subtracted out
        ROIs             % Nucleus mask (cell array of [x,y] polygon arrays)
        Parameters       % User-defined struct
        FileName         % List of loaded CSV files
        IsCulled  = false  % Flag: culling performed
        IsRelative = false % Flag: relative motion removed
        PixelSize     double = 0.165  % µm/px; camera constant; overridable via TOML
        FrameInterval double = NaN    % seconds; experiment-specific; must be set via TOML
        ParamsFile    char   = ''     % path to TOML file if loaded via readParams()
        SourceType    char   = 'csv'  % 'csv' or 'smd'
    end

    methods
        %% Constructor
        function obj = TrajectoryWrapper(varargin)
            % Optional name-value initialization.
            %   'Parameters'    - struct of pipeline parameters
            %   'PixelSize'     - µm/px (default 0.165)
            %   'FrameInterval' - seconds between frames (default NaN = not set)
            p = inputParser;
            addParameter(p, 'Parameters',    struct());
            addParameter(p, 'PixelSize',     0.165);
            addParameter(p, 'FrameInterval', NaN);
            parse(p, varargin{:});

            obj.Parameters    = p.Results.Parameters;
            obj.PixelSize     = p.Results.PixelSize;
            obj.FrameInterval = p.Results.FrameInterval;
            obj.RawTracks     = {};
            obj.CulledTracks  = {};
            obj.RelativeTracks = {};
        end

        %% Load CSV files
        function obj = readData(obj, csvFile)
            % READDATA  Load trajectories from a CSV file.
            %
            % x,y coordinates from the Python tracker are in pixels and are
            % converted to µm on load using PixelSize.  Set PixelSize (via the
            % constructor or readParams) before calling readData.
            tbl  = readtable(csvFile, 'PreserveVariableNames', true);
            data = table2array(tbl(:, [1, 2, 16, 11, 18]));
            tracks = {};
            [~, idx] = sort(data(:,5));
            data = data(idx,:);
            trackIDs = unique(data(:,5));
            kk = 1;
            for itrack = 1:length(trackIDs)
                idx = find(data(:,5) == trackIDs(itrack));
                if length(idx) > obj.Parameters.minTrackLength
                    t        = data(idx, 1:4);
                    t(:,1:2) = t(:,1:2) * obj.PixelSize;   % px → µm
                    tracks{kk} = t;
                    kk = kk + 1;
                end
            end
            obj.RawTracks  = tracks;
            obj.FileName   = csvFile;
            obj.IsCulled   = false;
            obj.SourceType = 'csv';
            fprintf('Loaded %d trajectories.\n', length(tracks));
        end

        %% Read experiment parameters from a TOML file
        function obj = readParams(obj, tomlFile)
            % READPARAMS  Load experiment-specific parameters from a TOML file.
            %
            %   Recognized keys (case-insensitive field names after toml.map_to_struct):
            %     pixel_size      → PixelSize (µm/px)
            %     frame_interval  → FrameInterval (seconds)
            %   All other keys are merged into obj.Parameters.
            %
            %   Call before readData() or fromSMD() to set experiment params.

            config = toml.read(tomlFile);
            pars   = toml.map_to_struct(config);
            obj.ParamsFile = tomlFile;

            if isfield(pars, 'pixel_size')
                obj.PixelSize = pars.pixel_size;
                pars = rmfield(pars, 'pixel_size');
            end
            if isfield(pars, 'frame_interval')
                obj.FrameInterval = pars.frame_interval;
                pars = rmfield(pars, 'frame_interval');
            end

            % Merge remaining fields into Parameters
            fields = fieldnames(pars);
            for f = fields(:)'
                obj.Parameters.(f{1}) = pars.(f{1});
            end
        end

        %% Load nucleus mask from file (CSV path)
        function obj = loadNucleusMask(obj)
            % Finds the matching CellPose .mat mask file from FileName pattern.
            expr   = '_s(\d{3})_';
            tokens = regexp({obj.FileName}, expr, 'tokens');
            sxxx   = cellfun(@(c) c{1}{1}, tokens, 'UniformOutput', false);
            pattern = ['*s*', sxxx{1}, '_002_mask.mat'];
            maskFiles = dir(['../', pattern]);
            if isempty(maskFiles)
                pattern   = ['*s*', sxxx{1}, '_mask.mat'];
                maskFiles = dir(['../', pattern]);
            end
            if ~isempty(maskFiles)
                mask = load(['../', maskFiles.name]);
                bdry = bwboundaries(mask.mask);
                for iroi = 1:length(bdry)
                    x = bdry{iroi}(:,1) * obj.PixelSize;   % px → µm
                    y = bdry{iroi}(:,2) * obj.PixelSize;
                    bdry{iroi} = [x y];
                end
                obj.ROIs = bdry;
            end
        end

        %% Cull tracks using ROI membership + gap-consistency filter
        function culled_tracks = cull_tracks(obj)
            % Keeps only tracks that lie predominantly inside an ROI and
            % have no gaps in the first minLengthBeforeGap frames.
            if isempty(obj.RawTracks)
                error('No tracks loaded. Call readData() or fromSMD() first.');
            end
            if isempty(obj.ROIs)
                error('No nucleus mask loaded. Call loadNucleusMask() first.');
            end

            minLengthBeforeGap = obj.Parameters.minLengthBeforeGap;

            % Delegate to shared utility (frame numbers are in column 3 for TW)
            culled_tracks = TrackUtils.cullTracksCore( ...
                obj.RawTracks, obj.ROIs, 3, minLengthBeforeGap);

            % Shift frame numbers from 0-based to 1-based
            for ii = 1:length(culled_tracks)
                culled_tracks{ii}(:,3) = culled_tracks{ii}(:,3) + 1;
            end

            obj.CulledTracks = culled_tracks;
            obj.IsCulled     = true;
        end

        %% Remove relative (drift) motion
        function remove_relative_motion(obj)
            if isempty(obj.CulledTracks)
                error('No tracks loaded. Cull tracks first.');
            end
            data   = {};
            tracks = {};
            for itrack = 1:numel(obj.CulledTracks)
                data{itrack} = obj.CulledTracks{itrack}(:, [1:2, 3, 5]);
            end
            for iroi = 1:length(obj.ROIs)
                k  = 0;
                ix = [];
                for itrack = 1:numel(data)
                    if data{itrack}(1, end) == iroi
                        k      = k + 1;
                        ix(k)  = itrack;
                    end
                end
                if length(ix) > 10
                    try
                        [forest, ~, forestInfo] = maxSpanningForest(data(ix), 10, 2);
                        [relTracks, ~]          = relativeTracksFromForest(data(ix), forest);
                        if forestInfo.numComponents > 1
                            fprintf('  ROI %d: %d components\n', iroi, forestInfo.numComponents);
                        end
                        if numel(relTracks) > 5
                            tracks = [tracks, relTracks]; %#ok<AGROW>
                        end
                    catch
                        disp('cannot do correction');
                        continue;
                    end
                end
            end
            obj.RelativeTracks = tracks;
            obj.IsRelative     = true;
        end

        %% Getters
        function tracks = getRawTracks(obj)
            tracks = obj.RawTracks;
        end

        function tracks = getRelativeTracks(obj)
            if ~obj.IsRelative
                warning('Relative motion not yet removed. Returning culled tracks.');
                tracks = obj.CulledTracks;
            else
                tracks = obj.RelativeTracks;
            end
        end

        function tracks = getCulledTracks(obj)
            if ~obj.IsCulled
                warning('Tracks have not been culled yet. Returning raw tracks.');
                tracks = obj.RawTracks;
            else
                tracks = obj.CulledTracks;
            end
        end

        function n = getNumRawTracks(obj)
            n = length(obj.RawTracks);
        end

        function n = getNumCulledTracks(obj)
            n = length(obj.CulledTracks);
        end

        function mask = getNucleusMask(obj)
            mask = obj.ROIs;
        end

        function params = getParameters(obj)
            params = obj.Parameters;
        end

        function ps = getPixelSize(obj)
            ps = obj.PixelSize;
        end

        function fi = getFrameInterval(obj)
            fi = obj.FrameInterval;
        end

        function st = getSourceType(obj)
            st = obj.SourceType;
        end

        function obj = setPixelSize(obj, val)
            % SETPIXELSIZE  Override PixelSize (µm/px) after construction.
            obj.PixelSize = val;
        end

        function obj = setFrameInterval(obj, val)
            % SETFRAMEINTERVAL  Set FrameInterval (seconds) after construction.
            obj.FrameInterval = val;
        end

        %% Utility: quick stats
        function summary(obj)
            fprintf('\n=== TrajectoryWrapper Summary ===\n');
            fprintf('Source type      : %s\n', obj.SourceType);
            fprintf('Raw tracks       : %d\n', obj.getNumRawTracks());
            if obj.IsCulled
                fprintf('Culled tracks    : %d (%.1f%%)\n', ...
                    obj.getNumCulledTracks(), ...
                    100 * obj.getNumCulledTracks() / max(obj.getNumRawTracks(), 1));
            else
                fprintf('Culled tracks    : not yet performed\n');
            end
            fprintf('Nucleus mask     : %s\n', ...
                ifelse(isempty(obj.ROIs), 'not loaded', 'loaded'));
            fprintf('Pixel size       : %.4f µm/px\n', obj.PixelSize);
            if isnan(obj.FrameInterval)
                fprintf('Frame interval   : NOT SET (call readParams or set explicitly)\n');
            else
                fprintf('Frame interval   : %.4f s\n', obj.FrameInterval);
            end
            fprintf('Parameters stored: %d fields\n', length(fieldnames(obj.Parameters)));
            fprintf('==================================\n\n');
        end
    end

    methods (Static)
        function tw = fromSMD(smd, varargin)
        % FROMSMD  Create a TrajectoryWrapper from a processed SMD object.
        %
        %   tw = TrajectoryWrapper.fromSMD(smd)
        %   tw = TrajectoryWrapper.fromSMD(smd, 'PixelSize', 0.16)
        %
        %   The SMD object must have been processed through at minimum:
        %     smd.localize()
        %     smd.track()
        %     smd.get_roi()
        %     smd.cull_tracks()
        %
        %   Properties set automatically from SMD:
        %     PixelSize     = smd.pixelsize
        %     FrameInterval = 1 / smd.frame_rate
        %     ROIs          = smd.ROIs
        %     IsCulled      = true
        %     SourceType    = 'smd'
        %
        %   Track format is converted from SMD Nx10 (µm) to TW Nx5 (pixels)
        %   via TrajectoryAdapter.smdTracksToPixels.

            p = inputParser;
            addParameter(p, 'PixelSize', smd.pixelsize);
            parse(p, varargin{:});

            tw               = TrajectoryWrapper();
            tw.PixelSize     = p.Results.PixelSize;
            tw.FrameInterval = 1 / smd.frame_rate;
            tw.ROIs          = smd.ROIs;
            tw.SourceType    = 'smd';

            % Convert SMD tracks (Nx10, already in µm) → TW Nx5 format
            if ~isempty(smd.tracks)
                twFmt           = TrajectoryAdapter.smdTracksToMicrons(smd.tracks);
                tw.RawTracks    = twFmt;
                tw.CulledTracks = twFmt;  % SMD culling is in-place; no separate pre-cull raw
                tw.IsCulled     = true;
            else
                tw.RawTracks    = {};
                tw.CulledTracks = {};
                tw.IsCulled     = true;
            end

            % Relative tracks if available
            if ~isempty(smd.relative_tracks)
                tw.RelativeTracks = TrajectoryAdapter.smdRelativeToTW(smd.relative_tracks);
                tw.IsRelative     = true;
            end

            fprintf('fromSMD: %d culled tracks imported (PixelSize=%.4f µm/px, dt=%.4f s).\n', ...
                length(tw.CulledTracks), tw.PixelSize, tw.FrameInterval);
        end
    end
end

% -------------------------------------------------------------------------
% Local helper
% -------------------------------------------------------------------------
function str = ifelse(condition, trueStr, falseStr)
    if condition, str = trueStr; else, str = falseStr; end
end
