classdef TrackUtils
% TRACKUTILS  Shared static utilities for trajectory culling and gap filtering.
%
%   Used by both TrajectoryWrapper (CSV path) and SMD/cull_tracks (image path)
%   to avoid duplicated culling logic.

    methods (Static)

        function [culled, rawIndices, startFrames] = cullTracksCore(tracks, ROIs, frameCol, minLen)
        % CULLTRACKSCORE  ROI filter + gap-consistency filter.
        %
        %   [culled, rawIndices, startFrames] = TrackUtils.cullTracksCore(tracks, ROIs, frameCol, minLen)
        %
        %   tracks   - cell array of NxM matrices; cols 1:2 are x,y in µm
        %   ROIs     - cell array of [x, y] polygon vertex arrays (µm coords)
        %   frameCol - column index that holds the frame number
        %   minLen   - minimum consecutive frames before the first gap
        %              (minLengthBeforeGap parameter)
        %
        %   Returns:
        %     culled      - cell of surviving tracks with ROI_ID appended as last column
        %     rawIndices  - 1xK integer vector; rawIndices(k) is the index of culled{k}
        %                   in the input tracks cell array
        %     startFrames - 1xK vector; startFrames(k) is the frame number at which
        %                   culled{k} begins (after gap-based start trimming)

            ntracks    = length(tracks);
            culled     = {};
            rawIdxPre  = zeros(1, ntracks);
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
                    culled{k}           = tracks{ii};
                    culled{k}(:,end+1)  = ROI_belong;
                    rawIdxPre(k)        = ii;
                    k = k + 1;
                end
            end

            rawIdxPre = rawIdxPre(1:k-1);

            [culled, keepMask, startFrames] = ...
                TrackUtils.filterFragmentsWithInitialGaps(culled, frameCol, minLen);
            rawIndices = rawIdxPre(keepMask);
        end

        % ------------------------------------------------------------------

        function [frags, keepMask, startFrames] = filterFragmentsWithInitialGaps(frags, frameCol, M)
        % FILTERFRAGMENTSWITHINITIALGAPS  Discard fragments with early gaps.
        %
        %   [frags, keepMask, startFrames] = TrackUtils.filterFragmentsWithInitialGaps(frags, frameCol, M)
        %
        %   frags      - cell array of NxK matrices
        %   frameCol   - column index holding the frame number
        %   M          - required number of consecutive frames at the start
        %
        %   Returns:
        %     frags       - surviving fragments, each trimmed to the valid start
        %     keepMask    - 1xN logical; true for each input fragment that survived
        %     startFrames - 1xS vector (S = nnz(keepMask)); frame number at the
        %                   start of each surviving fragment after trimming

            n           = length(frags);
            validFrags  = {};
            keepMask    = false(1, n);
            startFrames = zeros(1, n);   % over-allocated; trimmed at end

            for i = 1:n
                frag = frags{i};
                if size(frag, 1) < M
                    continue;
                end

                frag   = sortrows(frag, frameCol);
                frames = frag(:, frameCol);

                startIdx = TrackUtils.findValidStart(frames, M);
                if ~isempty(startIdx)
                    validFrags{end+1}  = frag(startIdx:end, :); %#ok<AGROW>
                    keepMask(i)        = true;
                    startFrames(i)     = frames(startIdx);
                end
            end

            frags       = validFrags;
            startFrames = startFrames(keepMask);
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
