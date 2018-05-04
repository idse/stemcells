clear all;
close all;

[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(testScriptPath);
addpath(genpath(testScriptPath));
addpath(genpath('/Users/idse/repos/mechanicalFeedback/simulation/general_functions'));
addpath(genpath('/Users/idse/repos/epitheliumAnalysis'));
addpath(genpath('/Users/idse/repos/keep_external/toolbox_graph'));
%compile

%%
%-----------------------------------------
% GENERATE INITIAL LATTICE AS DISK
%-----------------------------------------

%addpath('./utils');

R = 500;
[X, Y] = meshgrid(-R:R,-R:R);
mask = X.^2 + Y.^2 < R.^2;

a = 14;
noise = 0.7*a;

[vertices, cc] = randomVoronoiLattice2(a, noise, mask);
g = GLattConversion2(cc,vertices,false);

clf
LatticePresentation2(g)

% create a CellLayer object called cellLayer
cellStateVars = {'C'};
bondStateVars = {};%'T'
nTimePts = 1;
cellLayer = CellLayer(nTimePts, cellStateVars, bondStateVars);

% initialize timepoint from Lattmin
t = 1;
cellLayer.initTime(t, 'Lattmin', g);

%%

cells = {cellLayer.cells{t}.bondInd};
bonds = cat(2,cat(1,cellLayer.bonds{t}.vertInd));%, cat(1,this.bonds{t}.cellInd));
verts = cellLayer.vertices{t};

nCells = numel(cells);
maxNN = 14;
faces = nan([nCells maxNN]);
CellCM = zeros([length(cells) 3]);
ColonyCM = mean(verts);

for ci = 1:length(cells)
    
    C = cells{ci};
    faces(ci,1:numel(C)) = bonds(C,1);
    
    CellCM(ci,:) = mean(verts(bonds(cells{ci},1),:)) - ColonyCM;
end

CellR = sqrt(sum(CellCM.^2,2));
maxR = max(CellR);

%%
colorIdx = zeros(size(CellR));
R1 = 0.75;
R2 = 0.65;
R3 = 0.5;
colorIdx(CellR > R1*maxR) = 1;
colorIdx(CellR > R2*maxR & CellR < R1*maxR) = 2;
colorIdx(CellR > R3*maxR & CellR < R2*maxR) = 3;
colorIdx(CellR < R3*maxR) = 4;

cmap = [0 1 0.3;
        1 0.9 0;
        1 0 0;
        0.2 0.2 1];
        
faceCol = cmap(colorIdx,:);

clf
patch('Faces',faces,'Vertices',verts,'FaceVertexCData',faceCol);
shading faceted
axis equal
axis off
drawnow

%% Bra/Sox2 picture only

CellR = sqrt(sum(CellCM.^2,2));
maxR = max(CellR);
colorIdx = zeros(size(CellR));
R1 = 0.75;
R2 = 0.65;
colorIdx(CellR > R1*maxR) = 1;
colorIdx(CellR < R1*maxR) = 2;

cmap = [1 0 0;
        0.2 0.2 1];
        
faceCol = cmap(colorIdx,:);

clf
patch('Faces',faces,'Vertices',verts,'FaceVertexCData',faceCol);
shading faceted
axis equal
axis off
drawnow
hold on
scatter(500,500,100,[0 0.7 0],'fill')
viscircles([500 500],maxR,'EdgeColor',[0 0.7 0],'LineWidth',3)
hold off

%% on sphere

% create spherical mesh from a pointcloud
R = 30;
N = 5000;
[X,Y,Z] = randSphere(N,R);
scatter3(X,Y,Z)
axis equal

TRI = delaunay(X,Y,Z);
TR = triangulation(TRI,[X Y Z]);
[B,V] = freeBoundary(TR);

% trisurf(B,V(:,1),V(:,2),V(:,3),0*V(:,1))
% axis equal
% axis off
% %shading interp

% % the triangles do seem to be properly oriented
% triNorm = cross(V(B(:,1),:) - V(B(:,3),:),V(B(:,2),:) - V(B(:,3),:));
% triNorm = bsxfun(@rdivide, triNorm, sqrt(sum(triNorm.^2,2)));
% R = V(B(:,3),:);
% R = bsxfun(@rdivide, R, sqrt(sum(R.^2,2)));
% sum(triNorm.*R,2);

% dual mesh
g = compute_dual_vertex_model(B,V);

% create a CellLayer object called cellLayer
cellStateVars = {'C'};
bondStateVars = {};%'T'
nTimePts = 1;
cellLayer = CellLayer(nTimePts, cellStateVars, bondStateVars);

% initialize timepoint from Lattmin
t = 1;
cellLayer.initTime(t, 'Lattmin', g);

% set state of cell 1 to true, rest to zero 
B0 = zeros([1 numel(cellLayer.cells{t})]);
B0(1)=1;
cellLayer.setCellState(t, B0');

%%

cells = {cellLayer.cells{t}.bondInd};
bonds = cat(2,cat(1,cellLayer.bonds{t}.vertInd));%, cat(1,this.bonds{t}.cellInd));
verts = cellLayer.vertices{t};

nCells = numel(cells);
maxNN = 14;
faces = nan([nCells maxNN]);
CellCM = zeros([length(cells) 3]);
ColonyCM = mean(verts);

for ci = 1:length(cells)
    
    C = cells{ci};
    faces(ci,1:numel(C)) = bonds(C,1);
    
    CellCM(ci,:) = mean(verts(bonds(cells{ci},1),:)) - ColonyCM;
end

CellR = sqrt(sum(CellCM.^2,2));
maxR = max(CellR);

%%
% visualize
maxZ = abs(max(CellCM(:,3)));
Bra = CellCM(:,3) > maxZ/2;

nCells = numel(g.cells);
maxNN = 14;
faces = nan([nCells maxNN]);
for ci = 1:length(g.cells) 
    C = g.cells{ci};
    faces(ci,1:numel(C)) = g.bonds(C,1);
end

fv = struct();
fv.Vertices = g.verts;
fv.Faces = faces;
fv.FaceVertexCData = zeros([size(Bra,1) 3]);
fv.FaceVertexCData(Bra,1) = 1;
fv.FaceVertexCData(~Bra,1) = 0.2;
fv.FaceVertexCData(~Bra,2) = 0.2;
fv.FaceVertexCData(~Bra,2) = 0.8;
%fv.FaceVertexCData = cellfun(@numel,g.cells');%0*fv.Vertices + 0.5;
fv.LineWidth = 1;

clf
patch(fv)
shading faceted
axis tight
axis equal
axis off
colormap hsv

hold on
%plot3([0 0], [0 0], 1.5*[-maxZ maxZ],'Color',[0 0.7 0],'LineWidth',3);
hold off

view([1 1 1])
camroll(30)

%set(gcf,'Color','w');

%camlight('left','local')
%lighting gouraud
