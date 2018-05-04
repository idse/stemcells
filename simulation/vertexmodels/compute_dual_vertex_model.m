function g = compute_dual_vertex_model(face,vertex)

% compute_dual_graph - compute the dual graph of a given triangulation
%
%   [A,vertex1] = compute_dual_graph(face,vertex);
%
%   'A' is the adjacency matrix of the abstract dual graph
%   (recall that this graph link togeter adjacent faces
%   in the triangulation).
%
%   'vertex' is optional, and if given, the position of the vertex
%   of the dual graph (contained in 'vertex1') will 
%   the centroids of the face's vertex positions.
%   
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    vertex = [];
end

bdry = unique(freeBoundary(triangulation(face,vertex)));

[vertex,face] = check_face_vertex(vertex,face);

nface = size(face,2);

% compute the center of the faces 
if nargin>1
    dverts = [   ...
    sum(reshape(vertex(1,face),[3 nface]), 1)/3; ...
    sum(reshape(vertex(2,face),[3 nface]), 1)/3; ...
    sum(reshape(vertex(3,face),[3 nface]), 1)/3 ];
else
    dverts = [];
end

% dual cells
ndbonds = 3*nface;
dbonds = zeros([ndbonds 4]);
fring = compute_face_ring(face);

dcells = cell([size(vertex,2) 1]);
dcellmask = false([size(vertex,2) 1]);

% for each face in the primal triangulation
for i=1:nface 
    
    % get the neighboring faces
    ring = fring{i};
    
    for j=1:length(ring)
            
        % determine primal edge c being intersected
        [~,ia] = setdiff(face(:,i), face(:,ring(j)));
        c = face(:,i);
        c(ia,:) = [];
        if ia == 2
            c = c([2 1]);
        end

        bi = 3*(i-1) + j;
        dbonds(bi,:) = [i, ring(j), c(2), c(1)];

        if c(2) > numel(dcells)
            dcells{c(2)} = bi;
        else
            dcells{c(2)} = [dcells{c(2)}, bi];
        end
    end
end


dcellmask(bdry) = true;
dcells(dcellmask) = [];

cidxmapping = zeros([size(vertex,2) 1]);
cidxmapping(~dcellmask) = 1:sum(~dcellmask);

dbonds(dbonds(:,3)~=0,3) = cidxmapping(dbonds(dbonds(:,3)~=0,3));
dbonds(dbonds(:,4)~=0,4) = cidxmapping(dbonds(dbonds(:,4)~=0,4));

% order the bonds
for  ci = 1:numel(dcells)
    if ~isempty(dcells{ci})
        X = cat(1, dbonds(dcells{ci},1:2));
        P = [1 0 0];
        for i = 2:numel(dcells{ci})
            P(i) = find(X(:,1) == X(P(i-1),2));
        end
        dcells{ci} = dcells{ci}(P);
    end
end

g = struct('cells', {dcells'},'bonds',dbonds,'verts',dverts','bdryBonds',[]);
% REMOVE INCOMPLETE CELLS CORRESPONDING TO TRI BOUNDARY VERTICES


