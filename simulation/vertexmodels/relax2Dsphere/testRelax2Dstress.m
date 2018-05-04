clear all;
close all;

[testScriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(testScriptPath);
addpath(genpath(testScriptPath));
addpath(genpath('/Users/idse/repos/mechanicalFeedback/simulation/general_functions'));
addpath(genpath('/Users/idse/repos/epitheliumAnalysis'));
%compile

%%
%-----------------------------------------
% GENERATE INITIAL LATTICE
%-----------------------------------------

%addpath('./utils');

N = 1;
a = 1;
noise = 0;

initialScale = 1;%1.41102

[vertices, cc] = randomVoronoiLattice(N, a, noise);
vertices = [initialScale*vertices 0*vertices(:,1)];

% cc = {[1 3 2 4]};
% vertices = [1 1 0; 1 -1 0; -1 -1 0; -1 1 0];

disp('Number of cells: ');
disp(length(cc));
if isempty(cc) 
    error('need at least one cell');
end

%([8 11 19 21]) ([8 21]) ([5 8 11 19 21 18 20])
g = GLattConversion2(cc,vertices,false);

% for area elasticity area of unit hexagonal lattice: sqrt(3)/2
a0 = (sqrt(3)/2);
areavar = 0;
g.A0 = a0*ones([numel(g.cells) 1]) + areavar*a0*randn([numel(g.cells) 1]); 
% % for inverse area term
% g.A0 = ones([numel(g.cells) 1]); 
g.kA0 = ones([numel(g.cells) 1]);% + rand([numel(g.cells) 1]);

% in the hexagonal lattice with unit lattice constant, the sides are
% 1/sqrt(3) => perimeter : 6 * 1/sqrt(3)
g.p0 = (6/sqrt(3))*ones([numel(g.cells) 1]); 
g.pT0 = ones([numel(g.cells) 1]);

g.l0 = (1/sqrt(3))*ones([size(g.bonds,1) 1]); 
% initialize with unit tensions
g.T0 = ones(size(g.bonds,1),1);
g.T = g.T0;

% clonal identities 
%----------------------
%g.clones = 1:numel(g.cells)';
g.clones = NaN*(1:numel(g.cells))';
%clsubset = [1:3 7:9 13:15];
%clsubset = [1:3 7:9 13:15];
%clsubset = 23;
%clsubset = [18 67 73 28 72 66 23];
clsubset = [];%[8 19];%[2];
g.clones(clsubset) = 1;%1:numel(clsubset);

actualArea = [];
for i = 1:numel(g.cells)
    X = g.verts(g.bonds(g.cells{i},1),1);
    Y = g.verts(g.bonds(g.cells{i},1),2);
    actualArea(i) = polyarea(X,Y);
end

% % set tension in 'clone' to zeros
% g.T0([g.cells{g.clones==1}]) = 0;
% g.T = g.T0;

% % set pressure to zero too
%g.A0(g.clones==1) = 4*g.A0(g.clones==1);
g.pT0(g.clones==1) = 0.5*g.pT0(g.clones==1);

clf
LatticePresentation2(g,struct('cellIndex',true))

% area, perimeter, bond tension, bond hookean
g.param = [1 0 0 1];

% this is clunky: 
ect = sumTension(g);    % for perimeter term
%ect = g.T;              % for hookean or string terms

%LatticePresentation2(gfinal, struct('transparent',true, 'edgeColor','b'))
LatticePresentation2(g, struct('transparent',true, 'edgeColor',ect, 'edgeColorAbsolute', true))

%%
%-----------------------------------------
% RELAX INITIAL LATTICE MECHANICALLY
%-----------------------------------------

% TODO: CREATE FIXED COLOR SCALE

%compile
clf
a = 1;
pT = 0.5; pH = 0; 
bT = 0; bH = 0; alpha = 0;%-1/2;

dispRange = [0.3 0.5];

g.param = [a pT pH bT bH alpha];
% tic
% gfinal = relaxEpithelium2D(g,1);
% toc 

%figure, LatticePresentation2(gfinal); %LatticePresentation(g,0,'b');
clf
LatticePresentation2(g,struct('cellIndex',true))

maxiter = 1000;
gfinal = relaxEpithelium2D(g, maxiter);
%mean(gfinal.T)

% this is clunky: 
ect = sumTension(gfinal);    % for perimeter term
%ct = gfinal.clones;
ct = gfinal.stress;
ct(isnan(ct))=0;

clf
%LatticePresentation2(gfinal, struct('transparent',true, 'edgeColor','b'))
LatticePresentation2(gfinal, struct('transparent', false, 'cellIndex',true, 'edgeColor',ect,...
 'edgeColorRange', dispRange, 'colorTable', ct, 'transparentColor', 0))

hold on  
%%
% %visualize forces
% s = 10;
% vidx = g.bonds([g.cells{[5 8 11 19 21 18 20]}],1);
%
% % % tension: red
% % quiver(gfinal.verts(:,1), gfinal.verts(:,2),...
% %                     s*gfinal.tForce(:,1), s*gfinal.tForce(:,2),0,'r');
% 
% % pressure: blue
% quiver(gfinal.verts(vidx,1), gfinal.verts(vidx,2),...
%                     s*gfinal.pForce(vidx,1), s*gfinal.pForce(vidx,2),0,'b');
%                 
% % perimeter magenta
% quiver(gfinal.verts(vidx,1), gfinal.verts(vidx,2),...
%                     s*gfinal.perimForce(vidx,1),s*gfinal.perimForce(vidx,2),0,'m');
% 


% SHOW THE FORCES ACTING ON A SINGLE CELL
%----------------------------------------------------
s= 10;
col = {'r','g','b'};
for i = 1%:3
    ci = [1 23 4];
    ci = ci(i);
Fc = gfinal.cellForce{ci};
quiver( gfinal.verts(gfinal.bonds(gfinal.cells{ci},1),1),...
        gfinal.verts(gfinal.bonds(gfinal.cells{ci},1),2),...
                     s*Fc(:,1),s*Fc(:,2),0,col{i},'LineWidth',2);
end
% 
% scatter(gfinal.verts(gfinal.bonds(gfinal.cells{1},1),1),...
%         gfinal.verts(gfinal.bonds(gfinal.cells{1},1),2));

% % net force: black
F = gfinal.tForce + gfinal.pForce + gfinal.perimForce;
quiver(gfinal.verts(:,1), gfinal.verts(:,2),...
                     s*F(:,1),s*F(:,2),0,'r');
             
V = gfinal.verts(gfinal.bonds(gfinal.cells{ci}(end),2),:);
scatter(V(1),V(2),'*g');

% WHAT IS WRONG WITH THE CELL FORCE ON THE STARRED VERTEX?

hold off
axis off

% how close is the net force to zero relative to the pressure?
relError = sqrt(sum(F.^2,2))./sqrt(sum(gfinal.pForce.^2,2));
%relError = relError(~isnan(relError));
max(relError)



% actualAreaFinal = [];
% for i = 1:numel(g.cells)
%     X = gfinal.verts(gfinal.bonds(gfinal.cells{i},1),1);
%     Y = gfinal.verts(gfinal.bonds(gfinal.cells{i},1),2);
%     actualAreaFinal(i) = polyarea(X,Y);
% end

% %% test hooks law was implemented correctly
% 
% T = sumTension(gfinal);
% l = sqrt(sum((gfinal.verts(gfinal.bonds(:,2),:) - gfinal.verts(gfinal.bonds(:,1),:)).^2, 2));
% Tm = (l - gfinal.l0);
% dT = T-Tm;
% 
% mean(dT) % avg error should be small

%% division

clf
gfinal = relaxEpithelium2D(g, maxiter);

ct(isnan(ct))=0;
options = struct('colorTable',ct);
LatticePresentation2(gfinal,options)

actualAreaFinal = [];
for i = 1:numel(gfinal.cells)
    X = gfinal.verts(gfinal.bonds(gfinal.cells{i},1),1);
    Y = gfinal.verts(gfinal.bonds(gfinal.cells{i},1),2);
    actualAreaFinal(i) = polyarea(X,Y);
end
sum(actualAreaFinal)

%%
clf
gfinal = relaxEpithelium2D(g, maxiter);
gfinal = divideCell(gfinal, 1, [1 0.1]);
gfinal = relaxEpithelium2D(gfinal, maxiter);
%mean(gfinal.T)
ct = gfinal.clones;
ct(isnan(ct))=0;
options = struct('colorTable',ct);
LatticePresentation2(gfinal,options)

actualAreaFinal = [];
for i = 1:numel(gfinal.cells)
    X = gfinal.verts(gfinal.bonds(gfinal.cells{i},1),1);
    Y = gfinal.verts(gfinal.bonds(gfinal.cells{i},1),2);
    actualAreaFinal(i) = polyarea(X,Y);
end
sum(actualAreaFinal)

s = 1;
ci = 1;
hold on
Fc = gfinal.cellForce{ci};
quiver( gfinal.verts(gfinal.bonds(gfinal.cells{ci},1),1),...
        gfinal.verts(gfinal.bonds(gfinal.cells{ci},1),2),...
                     s*Fc(:,1),s*Fc(:,2),0,col{i},'LineWidth',2);
hold off

hold on
s = 1000;
% % net force: black
F = gfinal.tForce + gfinal.pForce + gfinal.perimForce;
quiver(gfinal.verts(:,1), gfinal.verts(:,2),...
                     s*F(:,1),s*F(:,2),0,'r');

V = gfinal.verts(gfinal.bonds(gfinal.cells{ci}(end),2),:);
scatter(V(1),V(2),'*b');
hold off

%%
% this is clunky: 
ect = sumTension(gfinal);    % for perimeter term
ct = gfinal.clones;
ct(isnan(ct))=0;
dispRange = [min(ect) max(ect)];

LatticePresentation2(gfinal, struct('transparent', true, 'cellIndex',true, 'edgeColor',ect,...
 'edgeColorRange', dispRange, 'colorTable', ct, 'transparentColor', 0))

actualAreaFinal = [];
for i = 1:numel(gfinal.cells)
    X = gfinal.verts(gfinal.bonds(gfinal.cells{i},1),1);
    Y = gfinal.verts(gfinal.bonds(gfinal.cells{i},1),2);
    actualAreaFinal(i) = polyarea(X,Y);
end

%%
% implement the isogonal transformation based stress trace
%-----------------------------------------------------------

nCells = numel(gfinal.cells);
stress = zeros([nCells 1]);

for ci = 1:nCells
    
    vidx = gfinal.bonds(gfinal.cells{ci},1);
    vidx2 = gfinal.bonds(gfinal.cells{ci},2);
    verts = gfinal.verts(vidx,:);

    edgevec = gfinal.verts(vidx2,:) - gfinal.verts(vidx,:);

    bonds = gfinal.bonds;
    edgeBonds = bonds(bonds(:,4) == 0,:);
    dualEdgeBonds = edgeBonds(:,[2 1 4 3]);
    bonds = cat(1,bonds,dualEdgeBonds);

    outbondsPre = bonds(bonds(:,3) ~= ci & bonds(:,4) ~= ci,:);
    outbonds = zeros([size(gfinal.verts,1) 4]);
    outbonds(outbondsPre(:,1),:) = outbondsPre;
    outbonds = outbonds(vidx,:);
    bidx = outbonds(:,1) > 0;
    outbonds = outbonds(bidx,:);

    outvec = gfinal.verts(outbonds(:,2),:) - gfinal.verts(outbonds(:,1),:);% 
    % hold on
    % 
    % scatter(gfinal.verts(vidx(2:5),1),gfinal.verts(vidx(2:5),2),'or');
    % 
    % quiver(verts(:,1),verts(:,2),outvec(:,1),outvec(:,2),0,'g','LineWidth',2);
    % quiver(verts(:,1),verts(:,2),edgevec(:,1),edgevec(:,2),0,'r','LineWidth',2);
    % 
    % hold off
    % axis off

    % normalize vectors
    edgevec = bsxfun(@rdivide, edgevec, sqrt(sum(edgevec.^2,2)));
    outvec = bsxfun(@rdivide, outvec, sqrt(sum(outvec.^2,2)));
    edgevecshift = circshift(edgevec,[1 0]);

    % sines of angles
    sinip = cross(outvec, edgevec(bidx,:));
    sinim = cross(outvec, edgevecshift(bidx,:));
    sinip = sinip(:,3);
    sinim = sinim(:,3);

    factor = sinim./circshift(sinip, [-1 0]);

    dilation = bsxfun(@times, outvec, factor);

    sigma = sum(sum(gfinal.cellForce{ci}(bidx,:).*dilation(:,1:2), 2))/gfinal.A(ci);
    if isnan(sigma) 
        sigma = 0;
    end
    stress(ci) = sigma;
    
    isogonal = 0*verts;
    isogonal(bidx, :) = 0.1*dilation;
    isogonal = isogonal + verts;

%     fv = struct();
%     fv.Vertices = isogonal;
%     fv.Faces = 1:6;
%     patch(fv, 'FaceColor','none')
end

clf
LatticePresentation2(gfinal, struct('transparent', false, 'cellIndex',true, 'edgeColor',ect,...
 'edgeColorRange', dispRange, 'colorTable', stress))

%% test T1 transitions

bi = 22;

clf
presOpts = struct('transparent',false, 'edgeColor',ect, 'cellIndex', true);
LatticePresentation2(gfinal, presOpts)
%'colorTable', g.clones, 
x = gfinal.bonds(bi,:);
hold on
scatter(gfinal.verts(x(1:2),1), gfinal.verts(x(1:2),2),'*k')
hold off

%%
gT1 = T1transition(gfinal, bi);
gT1 = relaxEpithelium2D(gT1, maxiter);
clf
LatticePresentation2(gT1, presOpts)

% gT1 = T1transition(gT1, bi);
% gT1 = relaxEpithelium2D(gT1, maxiter);
% clf
% LatticePresentation2(gT1, presOpts)

% T1s don't quite square to the identity
% FIX? (could be irrelevant index permutation)
% DEFINITELY ALLOW T1s on BOUNDARY
all(gfinal.bonds == gT1.bonds);

hold on
ci = 17;
gx = gT1;
x = gx.verts(gx.bonds(gx.cells{ci},1),1);
y = gx.verts(gx.bonds(gx.cells{ci},1),2);
u = gx.verts(gx.bonds(gx.cells{ci},2),1) - x;
v = gx.verts(gx.bonds(gx.cells{ci},2),2) - y;
quiver(x,y,u,v,0,'g');
hold off

%% 
%-----------------------------------------
% RUN GROWTH SIMULATION
%-----------------------------------------

a = 2*10;
noise = 0.5*a;
R = 100;
[X,Y] = meshgrid(-R:R,-R:R);
mask = X.^2 + 1.5*Y.^2 < R.^2;
[vertices, cc] = randomVoronoiLattice2(a, noise, mask);
vertices = vertices/a;
vertices = bsxfun(@minus, vertices, mean(vertices));
ginit = GLattConversion2(cc,vertices,false);
clf
LatticePresentation2(ginit);

%%
%load('initialLatticeGrown65');
%ginit = gnew;
%ginit = g;

load('disklike1');
cloneindex = 77;

%load('disklikeSmall');
% cloneindex = 54;

% parameters
%-----------------------

% mechanics: 
% 1) area
% 2) perimeter tension, 3) perimeter hookean
% 5) bond tension 6) bond hookean
% 7) alpha

a = 1;
pT = 0.5; pH = 0; 
bT = 0; bH = 0;
alpha = 0;

param = [a pT pH bT bH alpha];

% growth rates: 
gamma = 0;
gammaClones = 1;
fbstrength = 0.05;

% dt ~> division rate
dt = 0.05; %0.1;

% numer of iterations in simulation
tmax = 150;

% initial data
%-----------------------

gnew = ginit;
%gnew = gfinal;
gnew.param = param;
gnew.A0 = ones([numel(gnew.cells) 1]);

% % temporary 
gnew.pT0    = ones([numel(gnew.cells) 1]);% + rand([numel(gnew.cells) 1]);
gnew.kA0    = ones([numel(gnew.cells) 1]);

% currently unused parameters
gnew.p0     = ones([numel(gnew.cells) 1]);
gnew.T0     = ones([size(gnew.bonds,1) 1]);
gnew.l0     = ones([size(gnew.bonds,1) 1]);

%gnew = gfinal;
gnew.clones = nan([numel(gnew.cells) 1]);
%gnew.clones([67 77 87 221 230 220 229]) = 1;
%gnew.clones(77) = 1;
gnew.clones(cloneindex) = 1;
% gnew.clones(34) = 1;
% gnew.clones(42) = 2;

maxiter = 10000;

gnew = relaxEpithelium2D(gnew,maxiter);

% %% test that tension is correctly assigned
% 
% disp('---')
% ci1 = 166; ci2 = 18;
% T1 = gnew.T(gnew.bonds(:,3) == ci1 & gnew.bonds(:,4) == ci2,:);
% T2 = gnew.T(gnew.bonds(:,4) == ci1 & gnew.bonds(:,3) == ci2,:);
% ect = sumTension(gnew);
% T1 + T2
% gnew.param(2)*(gnew.pT0(ci1) + gnew.pT0(ci2))
% ect(gnew.bonds(:,3) == ci1 & gnew.bonds(:,4) == ci2)

% tension display range
% Tmin = mean(gnew.T);
% Tmax = 3*mean(gnew.T);

Tmin = gnew.param(2)*2*mean(gnew.pT0) - 0.2;
Tmax = gnew.param(2)*2*mean(gnew.pT0) + 0.2;
TdispRange = [Tmin Tmax];
pDispRange = 2*[-Tmax Tmax];

% visualize inital frame
%-----------------------

% MESSED UP EDGE VISUALIZATION
% WHY DO IRRELEVANT NAN PARAMETERS MESS THINGS UP?

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4) scrsz(4)])

cTable = 0.3*ones([numel(gnew.cells) 1]);
cTable(1) = 0.5;

stress = computeDilationStress(gnew);

ect = sumTension(gnew); % edge color table
CT = gnew.clones; %stress;
CT(isnan(CT)) = 0;
options = struct(   'edgeColor',ect,...
                    'cellIndex',false,...
                    'transparent',false,...
                    'colorTable', CT,...
                    'edgeColorRange', TdispRange,...
                    'pDispRange', pDispRange);
%options = rmfield(options,'colorTable');
options.transparentColor = 0;
L = 20;
options.L = L;

clf
visualizeTension(gnew, options);

% legend
text(-L,-L,'t = 0');
info = ['area = ' num2str(a) '; perimHook = ' num2str(pH) ...
        '; perimTension = ' num2str(pT) ...
        '; bondHook = ' num2str(bH) '; bondTension = ' num2str(bT) ...
        '; \alpha = ' num2str(alpha)...
        '; \gamma_c/\gamma = ', num2str(gammaClones/gamma)];
    
text(-L,L,info);
text(-L,L - 1,['div rate = ' num2str(dt)]);
fname = fullfile(testScriptPath,'results','tension_frame0.tif');
saveas(gcf,fname);

%%
clf
visualizePressure(gnew, options);

% legend
text(-L,-L,'t = 0');
text(-L,L,info);
text(-L,L - 1,['div rate = ' num2str(dt)]);
fname = fullfile(testScriptPath,'results','pressure_frame0.tif');
saveas(gcf,fname);

% %visualize net force
% s = 1;
% hold on  
% someg = gnew;
% F = someg.tForce + someg.pForce;
% quiver(someg.verts(:,1), someg.verts(:,2),s*F(:,1),s*F(:,2),0,'k');
% hold off
    
% DELETE THIS
t = 1;
gtime = {};

%%
% run loop of division and relaxation
%--------------------------------------

%profile on
tmax = 250;
for t  = 1:tmax
    
    t, tic

    gcur = gnew;
    Ncells = numel(gcur.cells);

    % division
    %----------------------------------------
 
    % divide
    divParam = struct('A0', 1, 'k', 2);
    [divIdx, divAxis] = decideDivide(gnew, dt, 'area', divParam);
    
    divIdx = find(divIdx);
    for i = 1:numel(divIdx)
        ci = divIdx(i);
        gnew = divideCell(gnew, ci, divAxis(ci,:));
    end
    
    % before growth mechanical relaxation, this is just the old stress
    % before anything happened, but consistent with new number of cells
    oldStress = gnew.stress; 
    
    % mechanical relaxation
    %----------------------------------------
    
    gnew = relaxEpithelium2D(gnew, maxiter);

    % check that force balance still holds
    %---------------------------------------
    
    someg = gnew;
    F = someg.tForce + someg.pForce;
    relError = sqrt(sum(F.^2,2))./sqrt(sum(someg.pForce.^2,2));
    max(relError)
    mean(relError);
    
    disp('update parameters');
    %----------------------------------------
    
    outside = gnew.bonds(:,4)==0;
    %outsideCells = unique(gnew.bond(outside,3));
    outsideCells = false([numel(gnew.cells) 1]);
    outsideCells(gnew.bonds(outside,3)) = true;
    
%     % bond length and tension according to classical active T
%     g = gcur;
%     Lbond = sqrt(sum((g.verts(g.bonds(:,2),:) - g.verts(g.bonds(:,1),:)).^2,2));
%     g = gnew;
%     Lbondnew = sqrt(sum((g.verts(g.bonds(:,2),:) - g.verts(g.bonds(:,1),:)).^2,2));
%     dLbond = Lbondnew - Lbond; 
%     for ci = 1:Ncells
%         gnew.pT0(ci) = gnew.pT0(ci) + sum(dLbond(gnew.cells{ci}))/sum(Lbondnew(gnew.cells{ci}));
%     end

    Tmin = 0.1;
	gnew.pT0 = max(Tmin, gnew.pT0 - fbstrength.*(~outsideCells).*(gnew.stress - oldStress));
    if any(gnew.pT0 < 0)
        error('negative tension');
    end

    % UPDATE PARAMETERS PER CELL
    overgrowth = ~isnan(gnew.clones);
    
    % increase the area term as uniform growth
    gnew.A0 = gnew.A0 + ( gamma*(~overgrowth) + gammaClones*overgrowth )*dt;

    
    % visualization
    %----------------------------------------
    
    if any(abs(gnew.verts(:)) > L)
        L = 2*L;
        options.L = L;
    end

    clf
    visualizeTension(gnew, options);
    text(-L,-L,['t = ' num2str(t)]);
    drawnow
    
    fname = fullfile(testScriptPath,'results',['tension_frame' num2str(t) '.tif']);
    saveas(gcf,fname);
    
    clf
    visualizePressure(gnew, options);
    text(-L,-L,['t = ' num2str(t)]);
    drawnow
    
    fname = fullfile(testScriptPath,'results',['pressure_frame' num2str(t) '.tif']);
    saveas(gcf,fname);
    
    gtime{t} = gnew;
    
    t, toc
end 

%profile viewer

%%
t = 51;
g = gtime{t};
%g.T = circshift(g.T,[1 0]);
[min(g.stress) max(g.stress)] 
clf
options.cellIndex = false;
options.transparent = true;
options.transparentColor = NaN;
options.colorTable = g.stress;
visualizeTension(g, options);


%%
% 
% save(fullfile(testScriptPath,'initialLatticeGrown100tension'), 'gnew');
% ginit = gnew;

%% combine frames

tstamp = datestr(clock,'ddmmyy_HHMM');

type = {'tension', 'pressure'};

for i = 1:2

    fname = fullfile(testScriptPath,'results',['growth' tstamp '_' type{i} '.avi']);
    writerObj = VideoWriter(fname);
    writerObj.FrameRate = 2;
    open(writerObj);

    tcut = 65;

    for t = 0:tcut

        fname = fullfile(testScriptPath,'results',[type{i} '_frame' num2str(t) '.tif']);
        im = imread(fname);
        writeVideo(writerObj,im);
    end

    close(writerObj);
end

%%
% also save final frame without clonal markers

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(4) scrsz(4)])

options.cellIndex = false;
options.edgeColor = sumTension(gnew);  %,'colorTable',gnew.A);
%options = rmfield(options,'colorTabl e');
options.transparent = true;
dT = 0.1;
options.edgeColorRange = [1-dT 1+dT];

LatticePresentation2(gnew, options)

hold on
cloneCells = find(~isnan(gnew.clones));
cloneBonds = [gnew.cells{cloneCells}]';
cloneEdge = cloneBonds(isnan(gnew.clones(gnew.bonds(cloneBonds,4))));
V1 = gnew.verts(gnew.bonds(cloneEdge,1),:);
V2 = gnew.verts(gnew.bonds(cloneEdge,2),:);
for i = 1:size(V1,1)
    plot([V1(i,1) V2(i,1)],[V1(i,2) V2(i,2)],'k');
end
hold off   
    
axis([-L L -L L])
axis off
info = ['\gamma_c/\gamma = ', num2str(gammaClones/gamma), ';     t = ' num2str(t)];
text(-L,-L,info);

fname = fullfile(testScriptPath,'results',['growth' tstamp '_finalTension.tif']);
saveas(gcf,fname);

%% test areas of cell

gX = gnew;
ci = 16;
N = numel(gX.cells{ci});
CM = repmat(mean(gX.verts(gX.bonds(gX.cells{ci},1),:)), [N 1]);
A = gX.verts(gX.bonds(gX.cells{ci},1),:) - CM;
B = circshift(gX.verts(gX.bonds(gX.cells{ci},1),:),[-1 0]) - CM;
cross(A,B)
CM(1,:)

%% color by area

figure,
area = [];
for i = 1:numel(gnew.cells)
    X = gnew.verts(gnew.bonds(gnew.cells{i},1),1);
    Y = gnew.verts(gnew.bonds(gnew.cells{i},1),2);
    area(i) = polyarea(X,Y);
end
options = struct('edgeColor','k', 'transparent',false, 'colorTable',area);
LatticePresentation2(gnew, options)
    