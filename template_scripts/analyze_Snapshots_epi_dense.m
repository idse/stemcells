clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

% THIS IS THE FIXED 8-WELL
%dataDir = '/Volumes/IdseData/Quyen/160627-8well-Act-pluripotency';
% dataDir = '/Users/Idse/data_tmp/160627-8well-Act-pluripotency';
% btfname = fullfile(dataDir, 'w1.btf');
% vsifile = fullfile(dataDir,'w1.vsi');

dataDir = '/Users/Idse/data_tmp/160708-8well-RIvsnoRI-EPI-2-cropped';
vsifile = fullfile(dataDir,'w2_E6-CHIR-Act.vsi');
% 
% btfname = fullfile(dataDir, 'w6_E6-CHIR-Act.btf');
% vsifile = fullfile(dataDir,'w6_E6-CHIR-Act.vsi');

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

meta.channelLabel = {'DAPI','', '', ''};
%meta.channelLabel = {'Bra','H2B', 'Smad4', 'Sox17'};
nucChannel = 1;

save(fullfile(dataDir,'metaData'),'meta');

%% load data for parameter check

n = 2;
m = 2;
ymin = n*2048;
xmin = n*2048;
ymax = ymin + m*2048;
xmax = xmin + m*2048;
preview = {};

% 2-3 times slower but saves us the btf trouble
series = 1;
tic
img_bf = bfopen_mod(vsifile,xmin,ymin,xmax-xmin+1,ymax-ymin+1,series);
for ci = 1:meta.nChannels
    preview{ci} = img_bf{1}{ci,1};
end
toc

%% check parameters for generating nuclear mask

% 3.5 microns is a good scale for the filters
opts.s = round(3.5/meta.xres);
opts.areacutoff = round(opts.s^2);
opts.normalizedthresh = 0.5;
opts.absolutethresh = 1000;

%im = preview{1}(2500:3000, 200:700);
im = preview{1};%(3000:end, 1:600);

profile on
nucseg = dirtyNuclearSegmentation(im, opts);
profile viewer

% visualize nuclear mask
newmaskedge = nucseg - imerode(nucseg,strel('disk',1));
impp = imadjust(mat2gray(im));
overlay = cat(3, impp, newmaskedge, impp);
imshow(overlay)

%% create full mask

im = bfopen_mod(vsifile,[],[],[],[],series, nucChannel);
im = im{1}{nucChannel,1};

tic
nucseg = dirtyNuclearSegmentation(im, opts);
toc

[~,barefname,~] = fileparts(vsifile);
imwrite(nucseg, fullfile(dataDir, [barefname '_nucseg.tif']), 'Compression', 'none');

%% process data

% options = struct(   'nuclearSegmentation',  nucseg,...
%                     'imopenBGsub',          50);
profile on
options = struct('nuclearSegmentation',  nucseg);
p = Position(meta.nChannels, vsifile, 1);
p.extractData(dataDir,nucChannel,options)
profile viewer

% save
[~,barefname,~] = fileparts(vsifile);
save(fullfile(dataDir, [barefname '_positions.mat']),'p');

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

