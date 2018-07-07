clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

dataDir = '/Users/idse/data_tmp/0_George/07022018 IWP2 Quants';

filenrs = 4716;
colRadius = 350;  % in micron
DAPIChannel = 1;
 
%% preprocess

param = {   'colRadiiMicron', colRadius, 'colMargin', 50,...
            'cleanScale', 30, 'adjustmentFactor', 0.15,...
            'channelLabel',{'DAPI','ISL1','SOX9','PAX6'} };

for i = 1%:numel(filenrs)

    vsifile = fullfile(dataDir, sprintf('Process_%d.vsi',filenrs(i)));
    colonies = processVsi(vsifile, dataDir, param);
    save(fullfile(dataDir,['colonies_' num2str(filenrs(i))]), 'colonies');
end

%% make radial profiles

load(fullfile(dataDir,'metaData.mat'),'meta');

i = 1;
load(fullfile(dataDir,['colonies_' num2str(filenrs(i))]), 'colonies');

figure,
doubleNormalize = true;
[nucAvgAllNormalized, r] = plotAveragesNoSegmentation(meta,...
                            colRadius, DAPIChannel, colonies, doubleNormalize);
saveas(gcf,fullfile(dataDir,['Process_' num2str(filenrs(i)) '_radialProfiles.png']));

%% look at variability between colonies

plotInterColonyRadialVariability(meta, colonies)
saveas(gcf,fullfile(dataDir,['Process_' num2str(filenrs(i)) '_variability.png']));

%%
normalized = true;
plotInterColonyRadialVariability(meta, colonies, [], normalized)
saveas(gcf,fullfile(dataDir,['Process_' num2str(filenrs(i)) '_DAPInormalized_variability.png']));


%% load colonies and show

i = filenrs(1);
load(fullfile(dataDir,['colonies_' num2str(i) '.mat'])); 

coli = 1;
colDir = fullfile(dataDir,['colonies_' num2str(i)]);
DAPI = colonies(coli).loadImage(colDir, DAPIChannel);
imshow(DAPI,[]);
