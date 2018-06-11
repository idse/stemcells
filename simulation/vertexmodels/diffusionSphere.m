clear all;
close all;

[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(testScriptPath);
addpath(genpath(testScriptPath));
addpath(genpath('/Users/idse/repos/mechanicalFeedback/simulation/general_functions'));
addpath(genpath('/Users/idse/repos/mechanicalFeedback/simulation/relax2Dsphere'));
addpath(genpath('/Users/idse/repos/epitheliumAnalysis'));
addpath(genpath('/Users/idse/repos/keep_external/toolbox_graph'));
%compile

%%
clf

% create spherical mesh from a pointcloud
R = 30;
N = 1000;
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

%%
% visualize

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
fv.FaceVertexCData = cellfun(@numel,g.cells');%0*fv.Vertices + 0.5;
fv.LineWidth = 2;

clf
patch(fv)
shading faceted
axis tight
axis equal
axis off
colormap hsv
view([1 1 1])

%% relax mechanically

% DOES NOT WORK YET

cd('/Users/idse/repos/mechanicalFeedback/simulation/relax2Dsphere');
compile

a = 1;
pT = 0.5; pH = 0; 
bT = 0; bH = 0; alpha = 0;%-1/2;
g.pT0 = ones(size(g.cells));
g.A0 = 10*ones(size(g.cells));

% irrelevant: remove
g.p0 = 0*g.pT0;
g.l0 = 0*g.pT0;
g.T0 = 0*g.pT0;
g.kA0 = 0*g.pT0;

g.param = [a pT pH bT bH alpha];

gnew = relaxEpithelium2D(g, 100);

% visualize

nCells = numel(g.cells);
maxNN = 14;
faces = nan([nCells maxNN]);
for ci = 1:length(g.cells) 
    C = g.cells{ci};
    faces(ci,1:numel(C)) = g.bonds(C,1);
end

fv = struct();
fv.Vertices = gnew.verts;
fv.Faces = faces;
fv.FaceVertexCData = cellfun(@numel,g.cells');%0*fv.Vertices + 0.5;
fv.LineWidth = 2;

clf
patch(fv)
shading faceted
axis tight
axis equal
axis off
colormap hsv
view([1 1 1])

%%
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
% visualize 
clf    
states = cat(1,cellLayer.cells{1}.state);
options = struct('cellIndex', false, 'colorTable', B0);
options = struct('cellIndex', false, 'colorTable', cellfun(@numel,g.cells));
%options = struct();
cellLayer.visualize(1, options)
axis off

%%
% L is a laplacian matrix on the dual lattice, so Ncells X Ncells

addpath(genpath('/Users/idse/repos/keep_external/toolbox_graph/'));

verts = V;
faces = B;
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
mask(B0) = false;

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
saveResult = true;%false;

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

for tidx = 1:numel(Tlin)
    
%     B = Ylin(tidx,1:n);
%     N = Ylin(tidx,n+1:2*n);
%     L = Ylin(tidx,2*n+1:end);
    
    CT = Ylin(tidx,:);
    colorsIdx = round(255*min(CT - minVal, maxVal-minVal)/(maxVal-minVal) + 1);
    faceCol = cmap(colorsIdx,:);

    clf
    patch('Faces',faces,'Vertices',verts,'FaceVertexCData',faceCol,'LineWidth',2);
    shading faceted
    %camlight left
    %lighting phong
    axis equal
    axis off
    view([10 -180]);
    drawnow

    set(gcf,'color','w');
    
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

