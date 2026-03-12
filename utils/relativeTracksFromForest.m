function [relTracks, trackComponentMap] = relativeTracksFromForest(tracks, forest)
%RELATIVETRACKSFROMFOREST  Relative trajectories from spanning forest
%
%  [relTracks, trackComponentMap] = relativeTracksFromForest(tracks, forest)
%
%  Computes child-minus-parent relative displacements for every component
%  in the spanning forest. Output format is identical to relativeTracksFromMST.
%
%  Inputs:
%    tracks - cell array of Nx4 matrices [x, y, frame, ROI]
%    forest - struct array from maxSpanningForest
%
%  Outputs:
%    relTracks        - cell array of [dx, dy, frame, ROI] per edge
%    trackComponentMap - Nx1 vector mapping each track to its component (0=isolated)

N = numel(tracks);
trackComponentMap = zeros(N, 1);

if isempty(forest)
    relTracks = {};
    return;
end

% ------------------------------------------------------------
% 1. Pre-process tracks: sorted by frame
% ------------------------------------------------------------
trk = struct('frames',{}, 'x',{}, 'y',{}, 'ROI',{});
for k = 1:N
    M = tracks{k};
    [frame, ord] = sort(M(:,3));
    trk(k).frames = frame;
    trk(k).x = M(ord,1);
    trk(k).y = M(ord,2);
    trk(k).ROI = M(1,end);
end

% ------------------------------------------------------------
% 2. Process each component
% ------------------------------------------------------------
relTracks = {};

for c = 1:numel(forest)
    nodes = forest(c).nodeIndices;
    edgeTab = forest(c).edgeTable;
    root = forest(c).rootNode;

    % Mark component membership
    trackComponentMap(nodes) = c;

    if height(edgeTab) == 0
        continue;
    end

    % Build adjacency list (using original track indices)
    adj = containers.Map('KeyType','int32','ValueType','any');
    for ni = 1:numel(nodes)
        adj(nodes(ni)) = [];
    end
    for i = 1:height(edgeTab)
        u = edgeTab.Node1(i);
        v = edgeTab.Node2(i);
        adj(u) = [adj(u); v];
        adj(v) = [adj(v); u];
    end

    % BFS from root
    parent = containers.Map('KeyType','int32','ValueType','int32');
    visited = containers.Map('KeyType','int32','ValueType','logical');
    for ni = 1:numel(nodes)
        visited(nodes(ni)) = false;
    end

    queue = root;
    visited(root) = true;

    while ~isempty(queue)
        u = queue(1);
        queue(1) = [];
        neighbors = adj(u);
        for i = 1:numel(neighbors)
            v = neighbors(i);
            if ~visited(v)
                visited(v) = true;
                parent(v) = u;
                queue = [queue; v]; %#ok<AGROW>
            end
        end
    end

    % Compute relative tracks: child - parent at common frames
    for ni = 1:numel(nodes)
        child = nodes(ni);
        if child == root || ~isKey(parent, child)
            continue;
        end
        par = parent(child);

        common = intersect(trk(par).frames, trk(child).frames);
        if isempty(common)
            continue;
        end

        [~, idxP] = ismember(common, trk(par).frames);
        [~, idxC] = ismember(common, trk(child).frames);

        dx = trk(child).x(idxC) - trk(par).x(idxP);
        dy = trk(child).y(idxC) - trk(par).y(idxP);

        relTracks{end+1} = [dx, dy, common, repmat(trk(child).ROI, size(common))]; %#ok<AGROW>
    end
end

end
