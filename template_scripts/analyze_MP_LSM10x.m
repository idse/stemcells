clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

mainDataDir = '/Users/idse/data_tmp/0_George/180504_timed iwp2 exp/';
dataDir = mainDataDir;  % fullfile(mainDataDir, 'LSM10x');

fnameformat = 'Image%.4d_01.oib';
filename = sprintf(fnameformat, 1);
colRadius = 350;

meta = MetadataMicropattern(fullfile(dataDir,filename));
meta.colRadiiMicron = colRadius;
meta.colRadiiPixel = round(meta.colRadiiMicron/meta.xres);

filenrs = {1:3 6:8 [12 13 14] 18:20 24:26 30:32}; 
meta.channelLabel = {'DAPI','PAX6','ISL1','SOX9'};

% filenrs = {4:5 9:11 15:17 21:23 27:29 33:35};
% meta.channelLabel = {'DAPI','SIX1','PAX3','SNAI1'};

meta.conditions = strcat({'72','66','60','48','36','24'},'h IWP2');

save(fullfile(dataDir,'metaData.mat'),'meta');

DAPIChannel = 1;

%% process (save DAPI channel MIP for segmentation)

% although findColonies can find multiple colonies in a single image,
% for the LSM there is one per image, which this code below assumes

close all;

colonies(numel([filenrs{:}])) = Colony;

for coli = [filenrs{:}]
  
    % cleanScale in micron
    param = {'DAPIChannel',DAPIChannel, 'colID',coli, 'adjustmentFactor', 0.2};
    
    filename = sprintf(fnameformat, coli);
    colony = processOib(filename, dataDir, param,'type','SIP');
    colonies(coli) = colony;
    colonies(coli).setID(coli);
    % setID overwrites filename, which was the right thing for epi, not
    % here
    colonies(coli).filename = filename; 
end

save(fullfile(dataDir,'colonies'), 'colonies');

%% show plot of different conditions side by side

load(fullfile(dataDir,'colonies'), 'colonies');

% first normalize by DAPI, then scale all profiles from 0 to 1
doubleNormalize = true; 

for i = 1:numel(meta.conditions)    
    coloniesCombined{i} = colonies(filenrs{i});
end

figure('Position',[0 0 1600 300]);
plotMultipleAveragesNoSegmentation(meta, colRadius, DAPIChannel,...
                                coloniesCombined, meta.conditions, doubleNormalize)
saveas(gcf,fullfile(dataDir,'radialProfilesCombined.png'));                        


%% show plot of different conditions side by side w/o DAPI  normalize

doubleNormalize = true;
figure('Position',[0 0 1600 300]);
plotMultipleAveragesNoSegmentation(meta, colRadius, [],...
                                coloniesCombined, meta.conditions, doubleNormalize)
saveas(gcf,fullfile(dataDir,'radialProfilesCombinedNoDAPI.png'));                        

%% overlay conditions per stain

figure('Position',[0 0 1600 400]);
plotConditionOverlayNoSegmentation(meta, colRadius, [],...
                                coloniesCombined, meta.conditions, doubleNormalize)
saveas(gcf,fullfile(dataDir,'radialProfilesConditionOverlay.png'));

%% load colonies and show side by side

load(fullfile(dataDir,'metaData.mat'),'meta');
load(fullfile(dataDir,'colonies'), 'colonies');

% representative colonies for each condition
repcols = [7 4 6 6 7 7; 1 1 1 1 1 1; 12 12 12 12 12 12]'; 

ci = 1; % sox9

figure('Position',[0 0 1600 1000]),
m = 3;
n = numel(filenrs);

% set LUT from image of my choice, for SOX9 the 24h iwp is the brightest
col = colonies(filenrs{end}(1));
img = col.loadImage(colDir, ci);
img = max(img,[],3);
Ilim = round(stretchlim(img)*(2^16-1));

for i = 1:n
    for j = 1:m
   
        ii = (j-1)*n + i;

        col = colonies(filenrs{i}(m));
        colDir = dataDir;
        img = col.loadImage(colDir, ci);
        img = max(img,[],3);
        b= col.boundingBox;
        img = img(b(3):b(4),b(1):b(2));
        %bg = imopen(img, strel('disk',15));
        %img = img - bg;
%         if j==1 && i == 1
%             Ilim = round(stretchlim(img)*(2^16-1));
%         end
        subplot_tight(m,n,ii)
        imshow(img,Ilim);
        title(meta.conditions{i})
    end
end
saveas(gcf,fullfile(dataDir,['compareColonies_' meta.channelLabel{ci} '.png']));

