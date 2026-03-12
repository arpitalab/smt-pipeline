classdef TrackUtils
% TRACKUTILS  Shared static utilities for trajectory culling and gap filtering.
%
%   Used by both TrajectoryWrapper (CSV path) and SMD/cull_tracks (image path)
%   to avoid duplicated culling logic.

    methods (Static)

        function culled = cullTracksCore(tracks, ROIs, frameCol, minLen)
        % CULLTRACKSCORE  ROI filter + gap-consistency filter.
        %
        %   culled = TrackUtils.cullTracksCore(tracks, ROIs, frameCol, minLen)
        %
        %   tracks   - cell array of NxM matrices; cols 1:2 are x,y in pixels
        %   ROIs     - cell array of [x, y] polygon vertex arrays (pixel coords)
        %   frameCol - column index that holds the frame number
        %   minLen   - minimum consecutive frames before the first gap
        %              (minLengthBeforeGap parameter)
        %
        %   Returns culled: cell of tracks that (a) lie predominantly inside
        %   an ROI and (b) pass the gap-consistency criterion.
        %   ROI identity is appended as the last column of each track matrix.

            ntracks = length(tracks);
            culled  = {};
            k = 1;

            for ii = 1:ntracks
                tmp = tracks{ii}(:, 1:2);
                in  = zeros(size(tmp,1), length(ROIs));

                for jj = 1:length(ROIs)
                    in(:,jj) = inpolygon(tmp(:,1), tmp(:,2), ...
                                         ROIs{jj}(:,1), ROIs{jj}(:,2));
                end

                maximum    = 0;
                ROI_belong = 0;
                for jj = 1:length(ROIs)
                    if nnz(in(:,jj)) > 0.5 * size(tmp,1)
                        maximum    = max(nnz(in(:,jj)), maximum);
                        ROI_belong = jj;
                    end
                end

                if ROI_belong > 0
                    culled{k}             = tracks{ii};
                    culled{k}(:, end+1)   = ROI_belong;
                    k = k + 1;
                end
            end

            culled = TrackUtils.filterFragmentsWithInitialGaps(culled, frameCol, minLen);
        end

        % ------------------------------------------------------------------

        function frags = filterFragmentsWithInitialGaps(frags, frameCol, M)
        % FILTERFRAGMENTSWITHINITIALGAPS  Discard fragments with early gaps.
        %
        %   frags = TrackUtils.filterFragmentsWithInitialGaps(frags, frameCol, M)
        %
        %   frags    - cell array of NxK matrices
        %   frameCol - column index holding the frame number
        %   M        - required number of consecutive frames at the start

            validFragments = {};
            for i = 1:length(frags)
                frag = frags{i};
                if size(frag, 1) < M
                    continue;
                end

                frag   = sortrows(frag, frameCol);
                frames = frag(:, frameCol);

                startIdx = TrackUtils.findValidStart(frames, M);
                if ~isempty(startIdx)
                    validFragments{end+1} = frag(startIdx:end, :); %#ok<AGROW>
                end
            end
            frags = validFragments;
        end

        % ------------------------------------------------------------------

        function startIdx = findValidStart(frames, M)
        % FINDVALIDSTART  Earliest index where next M frames are consecutive.
        %
        %   startIdx = TrackUtils.findValidStart(frames, M)
        %
        %   frames - sorted column vector of frame numbers
        %   M      - required run length

            if length(frames) < M
                startIdx = [];
                return;
            end

            if all(diff(frames(1:M)) == 1)
                startIdx = 1;
                return;
            end

            diffs  = diff(frames);
            gapIdx = find(diffs > 1, 1, 'first');

            if isempty(gapIdx)
                if all(diff(frames(1:M)) == 1)
                    startIdx = 1;
                else
                    startIdx = [];
                end
            else
                startIdx = TrackUtils.findValidStart(frames(gapIdx+1:end), M);
                if ~isempty(startIdx)
                    startIdx = gapIdx + startIdx;
                end
            end
        end

    end
end
