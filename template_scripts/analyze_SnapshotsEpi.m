clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

% THIS IS THE FIXED 8-WELL
dataDir = '/Volumes/IdseData/160416-EndDIff-Sox17nS4cells-epi';

btfname = fullfile(dataDir, 'endo_w1.btf');
vsifile = fullfile(dataDir,'Image750.vsi');

% size limit of data chunks read in
maxMemoryGB = 4;

% TODO : local normalization by DAPI
% TODO : background subtraction by large scale Fourier transform or
% TODO : make objects save as struct to make sure we don't lose data       

%%
% metadata
%----------------

metaDataFile = fullfile(dataDir,'metaData.mat');

% read metadata
if exist(metaDataFile,'file')
    disp('loading previously stored metadata');
    load(metaDataFile);
elseif exist('vsifile','var') && exist(vsifile,'file')
    disp('extracting metadata from vsi file');
    meta = Metadata(vsifile);
    
% or pretty much enter it by hand
else
    meta = Metadata();
    h = Tiff(btfname);
    meta.ySize = h.getTag('ImageLength');
    meta.xSize = h.getTag('ImageWidth');
    h.close;
    meta.nChannels = 4;
    meta.channelNames = {'DAPI','GFP','RFP','CY5'};
    meta.xres = 0.650/2;
    meta.yres = meta.xres;
end

% manually entered metadata
%------------------------------

meta.channelLabel = {'Sox17ab','H2B', 'Sox17', 'eomes'};
%meta.channelLabel = {'Bra','H2B', 'Smad4', 'Sox17'};
nucChannel = 2;

save(fullfile(dataDir,'metaData'),'meta');

%% load data for parameter check

p = Position(meta.nChannels, btfname, 1);
sox17ab =  p.loadImage(dataDir, 1);
eomes =  p.loadImage(dataDir, 4);
sox17 =  p.loadImage(dataDir, 3);
seg = p.loadSegmentation(dataDir, nucChannel);
nuc = p.loadImage(dataDir, nucChannel);

%% check for good clean up parameters

idx = 2048:3*2048;
s17piece = sox17(idx,idx);
s17abpiece = sox17ab(idx,idx);
eomespiece = eomes(idx,idx);
nucpiece = nuc(idx,idx);
segpiece = seg(idx,idx);

cleanupOptions = struct('separateFused', true,...
                    'clearBorder',true, 'minAreaStd', 1, 'minSolidity',0,...
                    'openSize', 3,'minArea',50);
                
nucmask = nuclearCleanup(segpiece, cleanupOptions);
maskedge = imdilate(nucmask,strel('disk',1))-segpiece;

figure, imshow(cat(3,imadjust(mat2gray(nucpiece)),maskedge,imadjust(mat2gray(s17piece))),[])

%% process data

options = struct('cleanupOptions', cleanupOptions,'imopenBGsub',50);
p = Position(meta.nChannels, btfname, 1);
p.extractData(dataDir,nucChannel,options)

% save
[~,barefname,~] = fileparts(btfname);
save(fullfile(dataDir, [barefname '2.mat']),'p');

%% load processed data

[~,barefname,~] = fileparts(btfname);
load(fullfile(dataDir, [barefname '2.mat']));
p.cellData 

%% scatter cell locations on top of image

s17bg = imopen(s17piece,strel('disk',50));

h2b = imadjust(mat2gray(nucpiece));
s17 = imadjust(mat2gray(s17piece - s17bg));
s17ab = imadjust(mat2gray(s17abpiece  - imopen(s17abpiece,strel('disk',50))));
eom = imadjust(mat2gray(eomespiece - imopen(eomespiece,strel('disk',50))));

figure, imshow(cat(3,h2b,maskedge,s17),[])

%%

figure, imshow(cat(3,h2b,maskedge,s17),[])

cidx = p.cellData.XY(:,1) > idx(1) & p.cellData.XY(:,1) < idx(end)...
            & p.cellData.XY(:,2) > idx(1) & p.cellData.XY(:,2) < idx(end);
        
% color by normalize sox17-rfp
C = mat2gray(p.cellData.nucLevel(cidx,3)./p.cellData.nucLevel(cidx,2));

hold on
scatter(p.cellData.XY(cidx,1)-idx(1),p.cellData.XY(cidx,2)-idx(1),100,C,'.')
hold off

%% scatter plot of sox17 staining vs sox17-rfp

h2bLevels = p.cellData.nucLevel(:,nucChannel) - double(min(nuc(:)));
s17Levels = p.cellData.nucLevel(:,3) - double(min(sox17(:)));
s17abLevels = p.cellData.nucLevel(:,1) - double(min(sox17ab(:)));

scatter(s17abLevels./h2bLevels, s17Levels./h2bLevels);
axis([0 0.5 0 4])

%% single histogram

bright = s17Levels./h2bLevels > 0.15;
figure, hist(p.cellData.nucLevel(bright,3)./p.cellData.nucLevel(bright,2),50)

%%

figure, imshow(cat(3,s17ab,eom + 0*maskedge,s17),[])
%C = mat2gray(p.cellData.nucLevel(cidx,3)./p.cellData.nucLevel(cidx,2));
%bright = cidx & mat2gray(p.cellData.nucLevel(:,3)./p.cellData.nucLevel(:,2)) > 0.05;
%bright = cidx & mat2gray(p.cellData.nucLevel(:,3)) > 0.03;
%bright = cidx & p.cellData.nucLevel(:,3) - min(p.cellData.nucLevel(:,3)) > 100;
bright = cidx & s17Levels > mean(s17Levels);

hold on
scatter(p.cellData.XY(bright,1)-idx(1),p.cellData.XY(bright,2)-idx(1),100,'.c')
hold off

hold on
scatter(p.cellData.XY(cidx & ~bright,1)-idx(1),p.cellData.XY(cidx & ~bright,2)-idx(1),100,'.r')
hold off

