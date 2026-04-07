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
    %
    % Storage model
    % -------------
    % After cull_tracks(), every track in RawTracks is Nx5:
    %   cols 1-4 : x (µm), y (µm), frame (1-based), SNR
    %   col  5   : ROI_ID  (> 0 = passed culling; -1 = rejected)
    %
    % Tracks that pass culling are trimmed to their valid start frame.
    % Rejected tracks keep their original data with col 5 set to -1.
    % getCulledTracks() is a simple filter: tracks where col5 > 0.
    %
    % Before cull_tracks() is called, RawTracks entries are Nx4 (no col 5).

    properties (SetAccess = private)
        RawTracks        % Cell array of track matrices; Nx4 before culling, Nx5 after
        RelativeTracks   % Cell array: tracks with global cell motion subtracted out
        ROIs             % Nucleus mask (cell array of [x,y] polygon arrays)
        Parameters       % User-defined struct
        FileName         % Path to the loaded CSV file
        IsCulled  = false  % Flag: cull_tracks() has been run
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

            obj.Parameters     = p.Results.Parameters;
            obj.PixelSize      = p.Results.PixelSize;
            obj.FrameInterval  = p.Results.FrameInterval;
            obj.RawTracks      = {};
            obj.RelativeTracks = {};
        end

        %% Load CSV files
        function obj = readData(obj, csvFile)
            % READDATA  Load trajectories from a CSV file.
            %
            % x,y coordinates from the Python tracker are in pixels and are
            % converted to µm on load using PixelSize.  Set PixelSize (via the
            % constructor or readParams) before calling readData.
            % Frame numbers are 0-based in the CSV (Python convention) and are
            % shifted to 1-based here.
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
                    t(:,3)   = t(:,3) + 1;                  % 0-based → 1-based frames
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
            % Finds the matching CellPose .mat mask file for the loaded CSV.
            % The mask is expected to be in the parent directory of the CSV file.
            expr   = '_s(\d{3})_';
            tokens = regexp(obj.FileName, expr, 'tokens');
            if isempty(tokens)
                warning('TrajectoryWrapper:noSampleToken', ...
                    'Could not parse sample ID from filename: %s', obj.FileName);
                return;
            end
            sxxx = tokens{1}{1};

            % Search for mask files relative to the CSV file location
            [csvDir, ~, ~] = fileparts(obj.FileName);
            if isempty(csvDir)
                csvDir = pwd;
            end
            searchDir = fullfile(csvDir, '..');

            pattern   = ['*s*', sxxx, '_002_mask.mat'];
            maskFiles = dir(fullfile(searchDir, pattern));
            if isempty(maskFiles)
                pattern   = ['*s*', sxxx, '_mask.mat'];
                maskFiles = dir(fullfile(searchDir, pattern));
            end

            if ~isempty(maskFiles)
                maskPath = fullfile(maskFiles(1).folder, maskFiles(1).name);
                mask     = load(maskPath);
                bdry     = bwboundaries(mask.mask);
                for iroi = 1:length(bdry)
                    x = bdry{iroi}(:,1) * obj.PixelSize;   % px → µm
                    y = bdry{iroi}(:,2) * obj.PixelSize;
                    bdry{iroi} = [x y];
                end
                obj.ROIs = bdry;
            else
                warning('TrajectoryWrapper:maskNotFound', ...
                    'No mask file found for sample %s in %s', sxxx, searchDir);
            end
        end

        %% Cull tracks using ROI membership + gap-consistency filter
        function culled_tracks = cull_tracks(obj)
            % Assigns ROI membership to every track and updates RawTracks in-place.
            % After this call, every RawTracks{k} is Nx5:
            %   col 5 > 0  : track passed (value = ROI_ID); trimmed to valid start
            %   col 5 = -1 : track rejected; original data preserved
            if isempty(obj.RawTracks)
                error('No tracks loaded. Call readData() or fromSMD() first.');
            end
            if isempty(obj.ROIs)
                error('No nucleus mask loaded. Call loadNucleusMask() first.');
            end

            minLengthBeforeGap = obj.Parameters.minLengthBeforeGap;

            [rawCulled, rawIndices, startFrames] = TrackUtils.cullTracksCore( ...
                obj.RawTracks, obj.ROIs, 3, minLengthBeforeGap);

            % Build lookup: original index → (ROI_ID, start frame)
            nRaw    = numel(obj.RawTracks);
            roiLookup   = zeros(nRaw, 1);    % 0 = rejected
            frameLookup = zeros(nRaw, 1);
            for k = 1:numel(rawIndices)
                roiLookup(rawIndices(k))   = rawCulled{k}(1, end);   % ROI_ID from last col
                frameLookup(rawIndices(k)) = startFrames(k);
            end

            % Rewrite RawTracks with col 5 appended
            newTracks = cell(nRaw, 1);
            for ii = 1:nRaw
                t = obj.RawTracks{ii};
                if roiLookup(ii) > 0
                    % Passed: trim to valid start, append ROI_ID
                    t = sortrows(t, 3);
                    t = t(t(:,3) >= frameLookup(ii), :);
                    newTracks{ii} = [t, repmat(roiLookup(ii), size(t,1), 1)];
                else
                    % Rejected: keep full original track, append -1
                    newTracks{ii} = [t, repmat(-1, size(t,1), 1)];
                end
            end

            obj.RawTracks = newTracks;
            obj.IsCulled  = true;

            culled_tracks = obj.getCulledTracks();
        end

        %% Remove relative (drift) motion
        function remove_relative_motion(obj)
            if ~obj.IsCulled
                error('No culled tracks available. Call cull_tracks() first.');
            end
            culled = obj.getCulledTracks();
            data   = {};
            tracks = {};
            for itrack = 1:numel(culled)
                data{itrack} = culled{itrack}(:, [1:2, 3, 5]);
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

        function tracks = getCulledTracks(obj)
            % Returns only the tracks that passed culling (col5 > 0).
            if ~obj.IsCulled
                warning('Tracks have not been culled yet. Returning raw tracks.');
                tracks = obj.RawTracks;
                return;
            end
            mask   = cellfun(@(t) t(1,5) > 0, obj.RawTracks);
            tracks = obj.RawTracks(mask);
        end

        function tracks = getRelativeTracks(obj)
            if ~obj.IsRelative
                warning('Relative motion not yet removed. Returning culled tracks.');
                tracks = obj.getCulledTracks();
            else
                tracks = obj.RelativeTracks;
            end
        end

        function n = getNumRawTracks(obj)
            n = length(obj.RawTracks);
        end

        function n = getNumCulledTracks(obj)
            if ~obj.IsCulled
                n = 0;
                return;
            end
            n = sum(cellfun(@(t) t(1,5) > 0, obj.RawTracks));
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

        %% Serialization
        function s = saveobj(obj)
            % Strip RawTracks and RelativeTracks from the saved representation.
            % The live in-memory object is NOT modified — only the serialized
            % copy written to disk loses the track data.
            % Use tc.reloadTracks() after loading to rebuild from original CSVs.
            s.FileName      = obj.FileName;
            s.PixelSize     = obj.PixelSize;
            s.FrameInterval = obj.FrameInterval;
            s.Parameters    = obj.Parameters;
            s.ParamsFile    = obj.ParamsFile;
            s.SourceType    = obj.SourceType;
            s.ROIs          = obj.ROIs;
            s.IsCulled      = obj.IsCulled;
            s.IsRelative    = obj.IsRelative;
            % RawTracks and RelativeTracks deliberately omitted.
        end

        %% Utility: quick stats
        function summary(obj)
            fprintf('\n=== TrajectoryWrapper Summary ===\n');
            fprintf('Source type      : %s\n', obj.SourceType);
            if obj.IsCulled && isempty(obj.RawTracks)
                fprintf('Raw tracks       : (not loaded — call tc.reloadTracks())\n');
            else
                fprintf('Raw tracks       : %d\n', obj.getNumRawTracks());
            end
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
        %   smdTracksToMicrons already produces Nx5 (x,y,frame,SNR,ROI_ID),
        %   which matches the post-cull RawTracks format directly.

            p = inputParser;
            addParameter(p, 'PixelSize', smd.pixelsize);
            parse(p, varargin{:});

            tw               = TrajectoryWrapper();
            tw.PixelSize     = p.Results.PixelSize;
            tw.FrameInterval = 1 / smd.frame_rate;
            tw.ROIs          = smd.ROIs;
            tw.SourceType    = 'smd';

            if ~isempty(smd.tracks)
                tw.RawTracks = TrajectoryAdapter.smdTracksToMicrons(smd.tracks);
            else
                tw.RawTracks = {};
            end
            tw.IsCulled = true;

            if ~isempty(smd.relative_tracks)
                tw.RelativeTracks = TrajectoryAdapter.smdRelativeToTW(smd.relative_tracks);
                tw.IsRelative     = true;
            end

            fprintf('fromSMD: %d culled tracks imported (PixelSize=%.4f µm/px, dt=%.4f s).\n', ...
                numel(tw.RawTracks), tw.PixelSize, tw.FrameInterval);
        end

        function obj = loadobj(s)
            % Reconstruct a TrajectoryWrapper from the struct written by saveobj.
            % RawTracks and RelativeTracks will be empty; call tc.reloadTracks()
            % to rebuild them from the original CSV files.
            obj = TrajectoryWrapper();
            if isstruct(s)
                obj.FileName      = s.FileName;
                obj.PixelSize     = s.PixelSize;
                obj.FrameInterval = s.FrameInterval;
                obj.Parameters    = s.Parameters;
                obj.ParamsFile    = s.ParamsFile;
                obj.SourceType    = s.SourceType;
                obj.ROIs          = s.ROIs;
                obj.IsCulled      = s.IsCulled;
                obj.IsRelative    = s.IsRelative;
            else
                % Already a TrajectoryWrapper (no saveobj was involved)
                obj = s;
            end
        end
    end
end

% -------------------------------------------------------------------------
% Local helper
% -------------------------------------------------------------------------
function str = ifelse(condition, trueStr, falseStr)
    if condition, str = trueStr; else, str = falseStr; end
end
