clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

% THIS IS THE FIXED 8-WELL
%dataDir = '/Volumes/IdseData/Quyen/160627-8well-Act-pluripotency';
% dataDir = '/Users/Idse/data_tmp/160627-8well-Act-pluripotency';
% btfname = fullfile(dataDir, 'w1.btf');
% vsifile = fullfile(dataDir,'w1.vsi');

dataDir = '/Users/Idse/data_tmp/160708-8well-RIvsnoRI-EPI-2-cropped';
filelist = {'w1_E6.vsi','w2_E6-CHIR-Act.vsi', 'w6_E6-CHIR-Act.vsi', 'w8_mTeSR.vsi'};

vsifile = fullfile(dataDir,filelist{1});
[~,barefname,~] = fileparts(vsifile);

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

meta.channelLabel = {'DAPI','sox17', 'sox2', 'nanog'};
%meta.channelLabel = {'DAPI','nanog', 'sox2', 'oct4'};
%meta.channelLabel = {'Bra','H2B', 'Smad4', 'Sox17'};
nucChannel = 1;

save(fullfile(dataDir,'metaData'),'meta');

%% load data for parameter check

n = 1;
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

%% save preview to check in FIJI

filename = fullfile(dataDir, [barefname '_preview.tif']);
imwrite(preview{1}, filename);
for ci = 2:meta.nChannels
    imwrite(preview{ci}, filename, 'WriteMode','append');
end

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
save(fullfile(dataDir, [barefname '_positions.mat']),'p');

%% scatter cell locations on top of part of preview to check result

DAPI = mat2gray(preview{1});
XY = p.cellData.XY;

inpreview =     XY(:,1) < xmax & XY(:,1) > xmin...
                & XY(:,2) < ymax & XY(:,2) > ymin;
XY = XY(inpreview,:);

imshow(DAPI,[])
hold on
scatter(XY(:,1) - xmin,XY(:,2) - ymin, '.r');
hold off

%% load processed data for all files in the file list

allData = {};
nucLevelAll = [];

for i = 1:numel(filelist)

    vsifile = fullfile(dataDir,filelist{i});
    [~,barefname,~] = fileparts(vsifile);

    load(fullfile(dataDir, [barefname '_positions.mat']));
    allData{i} = p;
    
    N = p.cellData.nucLevel(:,nucChannel);
    nucLevel = bsxfun(@rdivide, p.cellData.nucLevel, N);
    nucLevelAll = cat(1, nucLevelAll, nucLevel);
end

%% determine limits for displaying intensities

tol = 0.001;
lim = {};

for i = 1:4
    A = nucLevelAll(:,i);
    lim{i} = stretchlim(mat2gray(A),tol)*(max(A) - min(A)) + min(A);
end

%% histogram of different channels
% for w2_E6-CHIR-Act.vsi we see low nanog levels but many saturating sox17 cells

for fi = 3%1:numel(filelist)
    
    vsifile = fullfile(dataDir,filelist{fi});
    [~,barefname,~] = fileparts(vsifile);
    
    p = allData{fi};    
    N = p.cellData.nucLevel(:,nucChannel);
    nucLevel = bsxfun(@rdivide, p.cellData.nucLevel, N);

    for channelIndex = 2%2:4

        bins = linspace(0,lim{channelIndex}(2),40);
        n = histc(nucLevel(:,channelIndex), bins);
        bar(bins,n)
        axis([bins(1) bins(end) 0 max(n)]);

    %    hist(p.cellData.nucLevel(:,channelIndex), 40)
        title(meta.channelLabel{channelIndex})
        filename = fullfile(dataDir, [barefname '_distribution_' meta.channelLabel{channelIndex}]);
        %saveas(gcf, filename);
        %saveas(gcf, [filename '.png']);
    end
end

%% look at image of DAPI, sox17, nanog
% for w2_E6-CHIR-Act.vsi notice that none of the high nanog staining is real

imshow(cat(3, imadjust(mat2gray(preview{1})), 0.5*mat2gray(preview{2}), mat2gray(preview{4})),[])

%% scatter plots of markers against each other

for fi = 1:numel(filelist)

    vsifile = fullfile(dataDir,filelist{fi});
    [~,barefname,~] = fileparts(vsifile);

    p = allData{fi};

    for channelIdx1 = 2:4;
        for channelIdx2 = channelIdx1+1:4

            N = p.cellData.nucLevel(:,nucChannel);
            A = p.cellData.nucLevel(:,channelIdx1)./N;
            B = p.cellData.nucLevel(:,channelIdx2)./N;

            scatter(A, B,'.');
            xlabel(meta.channelLabel{channelIdx1})
            ylabel(meta.channelLabel{channelIdx2})
            title(filelist{fi},'Interpreter','none');
            
            axis([0 lim{channelIdx1}(2) 0 lim{channelIdx2}(2)]);

            filename = fullfile(dataDir, [barefname '_scatter_' meta.channelLabel{channelIdx1} '_' meta.channelLabel{channelIdx2}]);
            saveas(gcf, filename);
            saveas(gcf, [filename '.png']);
        end
    end
end

%% determine fraction of sox17+ etc

cutoff = 2;
channelIdx = 2;

for fi = 1:numel(filelist)
    
    p = allData{fi};
    
    N = p.cellData.nucLevel(:,nucChannel);
    A = p.cellData.nucLevel(:,channelIdx)./N;

    disp([filelist{fi} ' : ' num2str(100*sum(A > cutoff)/numel(A),2)]);
end
