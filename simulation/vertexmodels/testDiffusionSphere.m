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
clf

% create spherical mesh from a pointcloud
R = 30;
N = 200;
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

%% visualize

nCells = numel(dcells);
maxNN = 14;
faces = nan([nCells maxNN]);
for ci = 1:length(dcells) 
    C = dcells{ci};
    faces(ci,1:numel(C)) = dbonds(C,1);
end

fv = struct();
fv.Vertices = dverts';
fv.Faces = faces;
fv.FaceVertexCData = 0*dverts' + 0.5;
fv.LineWidth = 2;

clf
patch(fv)
shading faceted
axis tight
axis equal
axis off

%% test dual lattic generation in 2D

N = 2;
a = 1;
noise = 0;

[vertices, cc] = randomVoronoiLattice(N, a, noise);
vertices = 0.5*[vertices 0*vertices(:,1)];
g = GLattConversion2(cc,vertices,false);%([5 8 11 19 21 18 20])

clf
LatticePresentation2(g)

% create a CellLayer object called cellLayer
cellStateVars = {'C'};
bondStateVars = {};%'T'
nTimePts = 1;
cellLayer = CellLayer(nTimePts, cellStateVars, bondStateVars);
cellLayer.initTime(1, 'Lattmin', g);

dualCL = cellLayer.dualize;

% convert dual cellLayer to mesh for standard mesh operations

VofB = cat(1,dualCL.bonds{1}.vertInd);
BofC = cat(1,dualCL.cells{1}.bondInd);
VofB = VofB(:,1);

faces = VofB(BofC);
verts = dualCL.vertices{1};

V = verts;
B = faces;

%%

fv = struct();
fv.Vertices = verts;
fv.Faces = faces;
fv.FaceVertexCData = 0*verts + 1;

clf
patch(fv)
shading faceted
axis equal

for i = 1:size(faces,1)
    X = mean(verts(faces(i,:),1:2));
    text(X(1) + 0.02,X(2) + 0.02,num2str(i),'Color','r');
end

% dual mesh
[dverts,dbonds,dcells]  = compute_dual_vertex_model(faces,verts);

i = dbonds(:,1);
j = dbonds(:,2);

i = i(i~=0);
j = j(j~=0);

xy = dverts';
X = [ xy(i,1) xy(j,1) repmat(NaN,size(i))]';
Y = [ xy(i,2) xy(j,2) repmat(NaN,size(i))]';
Z = [ xy(i,3) xy(j,3) repmat(NaN,size(i))]';

hold on
plot3(X, Y, Z,'LineWidth',2,'Color','k')
axis equal

for i = 1:size(verts,1)
	text(verts(i,1) + 0.1,verts(i,2),num2str(i))
end

hold off
axis off

%% check dual lattice is constructed correctly by using patch

nCells = numel(dcells);
maxNN = 14;
dfaces = nan([nCells maxNN]);
for ci = 1:length(dcells) 
    C = dcells{ci};
    dfaces(ci,1:numel(C)) = dbonds(C,1);
end

fv = struct();
fv.Vertices = dverts';
fv.Faces = dfaces;
fv.FaceVertexCData = 0*dverts' + 1;
fv.LineWidth = 2;

clf
patch(fv)
shading faceted
axis equal
axis off

%%
%-----------------------------------------
% GENERATE INITIAL LATTICE AS DISK
%-----------------------------------------

%addpath('./utils');

R = 500;
[X, Y] = meshgrid(-R:R,-R:R);
mask = X.^2 + Y.^2 < R.^2;

a = 24; %14
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

% outside cells will be source of bmp
bci = cat(1,cellLayer.bonds{1}.cellInd);
outside = bci(:,2) == 0;
outsideCells = unique(bci(outside,1));

% set C 
B0 = zeros([1 numel(cellLayer.cells{t})]);
B0(outsideCells)=1;
cellLayer.setCellState(t, B0');

%%
% visualize 
clf    
states = cat(1,cellLayer.cells{1}.state);
options = struct('cellIndex', false, 'colorTable', states);
%options = struct();
cellLayer.visualize(1, options)

%%
%-----------------------------------------
% GENERATE INITIAL SIMPLE
%-----------------------------------------

N = 2;
a = 1;
noise = 0;

[vertices, cc] = randomVoronoiLattice(N, a, noise);
vertices = 0.5*[vertices 0*vertices(:,1)];
g = GLattConversion2(cc,vertices,false);%([5 8 11 19 21 18 20])

clf
LatticePresentation2(g)

% create a CellLayer object called cellLayer
cellStateVars = {'C'};
bondStateVars = {};%'T'
nTimePts = 1;
cellLayer = CellLayer(nTimePts, cellStateVars, bondStateVars);

% initialize timepoint from Lattmin
clf
t = 1;
cellLayer.initTime(t, 'Lattmin', g);
cellLayer.visualize(t);

% set center BMP 
B0 = zeros([1 numel(cellLayer.cells{t})]);
B0(2)=1;
cellLayer.setCellState(t, B0');

%%
% visualize 
clf    
states = cat(1,cellLayer.cells{1}.state);
options = struct('cellIndex', false, 'colorTable', states);
%options = struct();
cellLayer.visualize(1, options)

%%
% diffusion between cells: construct dual lattice first

dualCL = cellLayer.dualize;

% %%
% figure, 
% %cellLayer.visualize(t, options)
% dualCL.visualize(t, struct('transparent',true,'edgeColor','blue'));
% 
% %%
% verts = cellLayer.vertices{1};
% hold on 
% scatter(verts(:,1), verts(:,2));
% hold off
% 
% %%

% convert dual cellLayer to mesh for standard mesh operations

VofB = cat(1,dualCL.bonds{1}.vertInd);
BofC = cat(1,dualCL.cells{1}.bondInd);
VofB = VofB(:,1);

faces = VofB(BofC);
verts = dualCL.vertices{1};

fv = struct();
fv.Vertices = verts;
fv.Faces = faces;
fv.FaceVertexCData = 0*verts + 1;

clf
patch(fv)
shading faceted
axis equal

% %% dual dual
% 
% [A,vertex1] = compute_dual_graph(faces,verts);
% vertex1 = vertex1';
% 
% hold on 
% scatter(vertex1(:,1),vertex1(:,2))
% hold off
% 
% % 'A' is the adjacency matrix of the abstract dual graph
% %   (recall that this graph link togeter adjacent faces
% %   in the triangulation).
% 
% adj_list = adjmatrix2list(A);
% 
% adj_list{2}

%%
% L is a laplacian matrix on the dual lattice, so Ncells X Ncells

addpath(genpath('/Users/idse/repos/keep_external/toolbox_graph/'));

tic
options = struct('symmetrize',false, 'normalize', true);
W = compute_mesh_weight(verts,faces,'distance', options);
Lap = W - diag(sum(W,2));
toc

%Lap = 6*bsxfun(@rdivide, Lap, sum(W>0,2)');
%Lap = Lap/6;


%% solve 

tic

n = numel(cellLayer.cells{1});

% parameters: diffusion
Db = 1;
Dn = 10*Db;
Dl = 100*Dn;

% parameters : nodal lefty binding
knl = 1; 

% parameters : nodal/lefty induction by bmp
k  = 0.1;

% initial conditions
B0 = zeros([1 n]) + 1; %B0(outsideCells)=1;
N0 = zeros([1 n]);
L0 = zeros([1 n]);

% boundary conditions: hold value fixed in boundary cells
% the boundary cells should not be considered real cells, just a trick
mask = true(size(B0'));
mask(outsideCells) = false;

% differential equations
% BMP, Nodal, Lefty
% production, diffusion, degradation
dBdt    = @(B,N,L)        + Db*Lap*B;
dNdt    = @(B,N,L) k*B    + Dn*Lap*N    - knl*L.*N;
dLdt    = @(B,N,L) k*N    + Dl*Lap*L    - knl*L.*N;

odefun = @(t,y) cat(1,  mask.*dBdt(y(1:n),y(n+1:2*n),y(2*n+1:end)),...
                        mask.*dNdt(y(1:n),y(n+1:2*n),y(2*n+1:end)),...
                        mask.*dLdt(y(1:n),y(n+1:2*n),y(2*n+1:end))...
                    );

% differential equations
% just nodal to test
dNdt    = @(N) k*B0' + Dn*Lap*N;
odefun = @(t,y) mask.*dNdt(y);
  
tmax = R.^2/(4*Dn);

tspan   = [0 tmax];                    
y0      = N0';
[T,Y]   = ode23(odefun,tspan,y0);

% interpolate on linear timescale
dt      = 1;
Tlin    = (tspan(1):dt:tspan(end))';
Ylin    = interp1(T,Y,Tlin);

toc

%% visualize

barename = 'nodalProductionDiffusion';
saveResult = false;

cells = {cellLayer.cells{t}.bondInd};
bonds = cat(2,cat(1,cellLayer.bonds{t}.vertInd));%, cat(1,this.bonds{t}.cellInd));
verts = cellLayer.vertices{t};

nCells = numel(cells);
maxNN = 14;
faces = nan([nCells maxNN]);
for ci = 1:length(cells) 
    C = cells{ci};
    faces(ci,1:numel(C)) = bonds(C,1);
end

minVal = 0;
maxVal = max(Ylin(end,:));
cmap = jet(256);

for tidx = 200%numel(Tlin)
    
%     B = Ylin(tidx,1:n);
%     N = Ylin(tidx,n+1:2*n);
%     L = Ylin(tidx,2*n+1:end);
    
    CT = Ylin(tidx,:);
    colorsIdx = round(255*min(CT - minVal, maxVal-minVal)/(maxVal-minVal) + 1);
    faceCol = cmap(colorsIdx,:);

    clf
    patch('Faces',faces,'Vertices',verts,'FaceVertexCData',faceCol,'LineWidth',2);
    shading faceted
    axis equal
    axis off
    drawnow
    
    if saveResult
        if tidx == 1
            export_fig(gcf,[barename '.tif'],'-native');
        else
            export_fig(gcf,[barename '.tif'],'-native','-append');
        end
    end
end

%%

% plot result radially
vertsCM = mean(verts);
vertsRel = verts;
vertsRel(:,1) = vertsRel(:,1) - vertsCM(1);
vertsRel(:,2) = vertsRel(:,2) - vertsCM(2);
vertsR = sqrt(sum(vertsRel.^2,2));

edges = [0:40:R R];
[nhist,bin] = histc(vertsR, edges);
bincen = conv(edges,[1 1]/2,'valid');

clf
hold on
for tidx = 1:50:tmax
    
    %N = Ylin(tidx,n+1:2*n);
    N = Ylin(tidx,:);

    NradAvg = zeros(size(nhist));
    for i = 1:numel(nhist)
        NradAvg(i) = mean(N(bin == i));
    end

    plot(bincen, NradAvg(1:end-1))
end
r = bincen;
analyticSoln = k*(R.^2 - r.^2)/(4*Dn);

A = analyticSoln(1)./NradAvg(1);
plot(r, analyticSoln/A,'r')
hold off

% normalization issue with the Laplacian?

%%
% NEXT: define equations for cell state dynamics

% not efficient to do in for loop cell for cell, make function in cellLayer
% does defeat the object oriented approach a little

%%

% generate lattice on some random geometry using meshlab?
% basically all I need is a function to generate the dual mesh in meshlab

