function [subtracks, parentID] = saSPT_splitTracks(tracks, splitSize)
% SASPT_SPLITTRACKS  Gap-aware track splitting for state-array SPT.
%
%   [subtracks, parentID] = saSPT_splitTracks(tracks, splitSize)
%
%   Splits each track into contiguous subtracks of length up to splitSize,
%   discarding any segment shorter than 3 points.  Unlike the original
%   splitTracks.m (Heckert et al. 2021), gaps in the frame numbers are
%   respected: a track is first cut at every frame discontinuity (i.e.
%   diff(frame) ~= 1) and only the resulting gap-free runs are reshaped
%   into M-frame chunks.
%
%   Inputs
%     tracks      : cell array of [x, y, frame, ...] matrices (µm, frame)
%     splitSize   : maximum subtrack length M
%
%   Outputs
%     subtracks   : 1 x N cell array of [x, y] matrices, each gap-free,
%                   length 3..splitSize
%     parentID    : 1 x N int vector mapping each subtrack to its source
%                   track index in the input cell array

    if nargin < 2 || isempty(splitSize)
        splitSize = 7;
    end

    nT       = numel(tracks);
    subtracks = cell(1, 0);
    parentID  = zeros(1, 0);
    k         = 1;

    for itrack = 1:nT
        t = tracks{itrack};
        if size(t, 1) < 3
            continue;
        end

        if size(t, 2) >= 3
            f = t(:, 3);
            % indices where the next frame is not contiguous
            cuts = [0; find(diff(f) ~= 1); size(t, 1)];
        else
            cuts = [0; size(t, 1)];
        end

        for iseg = 1:numel(cuts)-1
            seg = t(cuts(iseg)+1 : cuts(iseg+1), 1:2);
            L   = size(seg, 1);
            if L < 3
                continue;
            end

            if L <= splitSize
                subtracks{k} = seg;
                parentID(k)  = itrack;
                k = k + 1;
            else
                numsplits = floor(L / splitSize);
                X1 = seg;
                dumx = reshape(X1(1:numsplits*splitSize, 1), splitSize, numsplits);
                dumy = reshape(X1(1:numsplits*splitSize, 2), splitSize, numsplits);
                for jj = 1:numsplits
                    subtracks{k} = [dumx(:, jj), dumy(:, jj)];
                    parentID(k)  = itrack;
                    k = k + 1;
                end
                tail = mod(L, splitSize);
                if tail > 2
                    subtracks{k} = X1(numsplits*splitSize+1 : end, :);
                    parentID(k)  = itrack;
                    k = k + 1;
                end
            end
        end
    end
end
