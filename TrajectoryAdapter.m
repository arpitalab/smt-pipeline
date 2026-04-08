classdef TrajectoryAdapter
% TRAJECTORYADAPTER  Static format-conversion utilities between SMD and
%                    TrajectoryWrapper track representations.
%
%   All tracks inside TrajectoryWrapper are stored in pixels.
%   These helpers translate the SMD Nx10 format and apply / reverse the
%   pixel-to-micron scaling that analysis methods expect.

    methods (Static)

        function twTracks = smdTracksToMicrons(smdTracks)
        % SMDTRACKSTOMICRONS  Reformat post-cull SMD tracks → TW Nx5 format (µm).
        %
        %   twTracks = TrajectoryAdapter.smdTracksToMicrons(smdTracks)
        %
        %   SMD post-cull track format (Nx10):
        %     col 1  : x position (µm)
        %     col 2  : y position (µm)
        %     col 3  : amplitude A
        %     col 4  : PSF width σ
        %     col 5  : frame number
        %     col 6  : SNR
        %     col 7  : iframe
        %     col 8  : neg log-likelihood
        %     col 9  : variance
        %     col 10 : ROI_ID
        %
        %   Output TW format (Nx5):
        %     col 1  : x position (µm)  — unchanged, already in µm
        %     col 2  : y position (µm)  — unchanged
        %     col 3  : frame number
        %     col 4  : SNR
        %     col 5  : ROI_ID

            twTracks = cell(size(smdTracks));
            for i = 1:length(smdTracks)
                t = smdTracks{i};
                if size(t, 2) >= 10
                    roi_id = t(:, 10);
                else
                    % get_roi() was not called — assign dummy ROI_ID of 1
                    roi_id = ones(size(t, 1), 1);
                end
                twTracks{i} = [t(:,1), ...   % x (µm)
                                t(:,2), ...   % y (µm)
                                t(:,5), ...   % frame
                                t(:,6), ...   % SNR
                                roi_id];      % ROI_ID (1 if no mask)
            end
        end

        % ------------------------------------------------------------------

        function twRelative = smdRelativeToTW(smdRelTracks)
        % SMDRELATIVETOTW  Validated passthrough: SMD relative tracks → TW format.
        %
        %   Both SMD and TW relative tracks share the format:
        %     [dx, dy, frame, ROI_ID]
        %
        %   This method asserts the column count matches and returns a copy.

            twRelative = cell(size(smdRelTracks));
            for i = 1:length(smdRelTracks)
                t = smdRelTracks{i};
                if size(t, 2) < 4
                    error('TrajectoryAdapter:smdRelativeToTW', ...
                          'Track %d has %d columns; expected at least 4 [dx,dy,frame,ROI].', ...
                          i, size(t,2));
                end
                twRelative{i} = t(:, 1:4);
            end
        end

        % ------------------------------------------------------------------

        function physTracks = pixelsToMicrons(twTracks, pixelSize)
        % PIXELSTOMICRONS  Scale x,y columns (1:2) from pixels to µm.
        %
        %   physTracks = TrajectoryAdapter.pixelsToMicrons(twTracks, pixelSize)
        %
        %   twTracks  - cell array of NxM track matrices (cols 1:2 are x,y in px)
        %   pixelSize - scalar µm/px conversion factor
        %
        %   Replaces the formerly hardcoded "* 0.165" in analysis methods.

            physTracks = cell(size(twTracks));
            for i = 1:length(twTracks)
                physTracks{i}       = twTracks{i};
                physTracks{i}(:,1)  = twTracks{i}(:,1) * pixelSize;
                physTracks{i}(:,2)  = twTracks{i}(:,2) * pixelSize;
            end
        end

    end
end
