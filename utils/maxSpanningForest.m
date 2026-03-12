function [forest, G, info] = maxSpanningForest(tracks, minOverlap, maxGap)
%MAXSPANNINGFOREST  Maximum spanning forest for disconnected track components
%
%  [forest, G, info] = maxSpanningForest(tracks, minOverlap, maxGap)
%
%  Like maxSpanningTreeFromTracksGapStrict, but handles multiple connected
%  components. Returns one MST per component so that ALL tracks participate
%  in drift correction, not just those reachable from node 1.
%
%  Inputs:
%    tracks     - cell array of Nx4 matrices [x, y, frame, ROI]
%    minOverlap - minimum real overlap frames for an edge (default 10)
%    maxGap     - max gap frames to fill for colocalization (default 2)
%
%  Outputs:
%    forest - struct array (one per component with >= 2 tracks)
%      .nodeIndices  - original track indices in this component
%      .edgeTable    - table(Edge, Node1, Node2, Weight) with original IDs
%      .rootNode     - track with most observed frames (BFS root)
%      .frameSpan    - [minFrame, maxFrame] across component
%      .numTracks    - number of tracks in component
%    G    - full overlap graph (all components)
%    info - struct with .numComponents, .isolatedTracks, .Freal

if nargin < 2, minOverlap = 10; end
if nargin < 3, maxGap = 2; end

N = numel(tracks);

% ------------------------------------------------------------
% 1. Pre-process
% ------------------------------------------------------------
trk = struct('frames',{}, 'x',{}, 'y',{}, 'ROI',{});
for k = 1:N
    M = tracks{k};
    [f,ord] = sort(M(:,3));
    trk(k).frames = f;
    trk(k).x = M(ord,1);
    trk(k).y = M(ord,2);
    trk(k).ROI = M(1,end);
end

% ------------------------------------------------------------
% 2. Pairwise real overlap
% ------------------------------------------------------------
Freal = zeros(N,N);

for i = 1:N-1
    fi = trk(i).frames;  xi = trk(i).x;  yi = trk(i).y;
    for j = i+1:N
        fj = trk(j).frames;  xj = trk(j).x;  yj = trk(j).y;

        framesUnion = union(fi, fj);
        if numel(framesUnion) < minOverlap, continue; end

        [xif, ~] = interpTrackWithGap(fi, xi, yi, framesUnion, maxGap);
        [xjf, ~] = interpTrackWithGap(fj, xj, yj, framesUnion, maxGap);

        % Real overlap: both originally present
        realMask_i = ismember(framesUnion, fi);
        realMask_j = ismember(framesUnion, fj);
        realCount = sum(realMask_i & realMask_j);

        if realCount >= minOverlap
            Freal(i,j) = realCount;
        end
    end
end
Freal = Freal + Freal.';

% ------------------------------------------------------------
% 3. Build graph
% ------------------------------------------------------------
[src, snk] = find(triu(Freal,1));
weights = Freal(sub2ind(size(Freal), src, snk));

if isempty(weights)
    G = graph([],[],[],N);
else
    G = graph(src, snk, weights, N);
end

% ------------------------------------------------------------
% 4. Find connected components
% ------------------------------------------------------------
compIDs = conncomp(G);
numComp = max(compIDs);

% Count observed frames per track (for root selection)
numObserved = zeros(N,1);
for k = 1:N
    numObserved(k) = numel(trk(k).frames);
end

% ------------------------------------------------------------
% 5. Build MST per component
% ------------------------------------------------------------
forest = struct('nodeIndices',{}, 'edgeTable',{}, 'rootNode',{}, ...
                'frameSpan',{}, 'numTracks',{});
isolatedTracks = [];
compIdx = 0;

for c = 1:numComp
    nodes = find(compIDs == c);
    if numel(nodes) < 2
        isolatedTracks = [isolatedTracks, nodes]; %#ok<AGROW>
        continue;
    end

    % Extract subgraph and compute MST
    Gsub = subgraph(G, nodes);
    GsubNeg = graph(Gsub.Edges.EndNodes(:,1), Gsub.Edges.EndNodes(:,2), ...
                    -Gsub.Edges.Weight, numel(nodes));
    Tsub = minspantree(GsubNeg);

    % Remap local node IDs back to original track indices
    localEdges = Tsub.Edges.EndNodes;
    origNode1 = nodes(localEdges(:,1))';
    origNode2 = nodes(localEdges(:,2))';
    origWeight = -Tsub.Edges.Weight;

    edgeTab = table((1:size(localEdges,1))', origNode1, origNode2, origWeight, ...
                    'VariableNames', {'Edge','Node1','Node2','Weight'});
    edgeTab = sortrows(edgeTab, 'Weight', 'descend');

    % Root = track with most observed frames in this component
    [~, bestLocal] = max(numObserved(nodes));
    rootNode = nodes(bestLocal);

    % Frame span
    allFrames = [];
    for ni = 1:numel(nodes)
        allFrames = [allFrames; trk(nodes(ni)).frames]; %#ok<AGROW>
    end
    frameSpan = [min(allFrames), max(allFrames)];

    compIdx = compIdx + 1;
    forest(compIdx).nodeIndices = nodes;
    forest(compIdx).edgeTable = edgeTab;
    forest(compIdx).rootNode = rootNode;
    forest(compIdx).frameSpan = frameSpan;
    forest(compIdx).numTracks = numel(nodes);
end

% ------------------------------------------------------------
% 6. Diagnostics
% ------------------------------------------------------------
info.numComponents = compIdx;
info.isolatedTracks = isolatedTracks;
info.Freal = Freal;

end

% ----------------------------------------------------------------
% Helper: fill gaps <= maxGap with last known position
% ----------------------------------------------------------------
function [xOut, yOut] = interpTrackWithGap(framesIn, xIn, yIn, framesOut, maxGap)
    xOut = nan(size(framesOut));
    yOut = nan(size(framesOut));

    [~, loc] = ismember(framesIn, framesOut);
    xOut(loc) = xIn;
    yOut(loc) = yIn;

    % Forward fill
    last = 0; gap = 0;
    for t = 1:numel(framesOut)
        if ~isnan(xOut(t))
            last = t; gap = 0;
        else
            gap = gap + 1;
            if gap <= maxGap && last > 0
                xOut(t) = xOut(last);
                yOut(t) = yOut(last);
            end
        end
    end

    % Backward fill start
    first = find(~isnan(xOut),1);
    if ~isempty(first)
        xOut(1:first-1) = xOut(first);
        yOut(1:first-1) = yOut(first);
    end
end
