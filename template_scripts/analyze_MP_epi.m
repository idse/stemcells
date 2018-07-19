clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

dataDir = '/Users/idse/data_tmp/0_George/07022018 IWP2 Quants';

filenrs = [4716 4719 4722 4725 4728 4731];
colRadius = 350;  % in micron
DAPIChannel = 1;
 
%% preprocess

param = {   'colRadiiMicron', colRadius, 'colMargin', 50,...
            'cleanScale', 30, 'adjustmentFactor', 0.15,...
            'channelLabel',{'DAPI','ISL1','SOX9','PAX6'},...
            'momentChannel',4, 'CMcutoff', 80};

for i = 1:numel(filenrs)

    vsifile = fullfile(dataDir, sprintf('Process_%d.vsi',filenrs(i)));
    colonies = processVsi(vsifile, dataDir, param);
    save(fullfile(dataDir,['colonies_' num2str(filenrs(i))]), 'colonies');
end

%% make radial profiles

load(fullfile(dataDir,'metaData.mat'),'meta');

for i = 1:numel(filenrs)
    load(fullfile(dataDir,['colonies_' num2str(filenrs(i))]), 'colonies');

    % make mean radial profiles
    figure,
    doubleNormalize = true;
    [nucAvgAllNormalized, r] = plotAveragesNoSegmentation(meta,...
                                colRadius, DAPIChannel, colonies, doubleNormalize);
    saveas(gcf,fullfile(dataDir,['Process_' num2str(filenrs(i)) '_radialProfiles.png']));
    
    % look at variability between colonies
    plotInterColonyRadialVariability(meta, colonies)
    saveas(gcf,fullfile(dataDir,['Process_' num2str(filenrs(i)) '_variability.png']));
    
    normalized = true;
    plotInterColonyRadialVariability(meta, colonies, [], normalized)
    saveas(gcf,fullfile(dataDir,['Process_' num2str(filenrs(i)) '_DAPInormalized_variability.png']));
end

%% combine means in single plot

conditions = strcat({'72','66','60','48','36','24'},'h IWP2');

for i = 1:numel(filenrs)    
    load(fullfile(dataDir,['colonies_' num2str(filenrs(i))]), 'colonies');
    coloniesCombined{i} = colonies;
end

figure('Position',[0 0 1600 300]),
plotMultipleAveragesNoSegmentation(meta, colRadius, DAPIChannel,...
                                coloniesCombined, conditions, doubleNormalize)
saveas(gcf,fullfile(dataDir,'radialProfilesCombined.png'));

%% load colonies and show

i = filenrs(1);
load(fullfile(dataDir,['colonies_' num2str(i) '.mat'])); 

coli = 1;
colDir = fullfile(dataDir,['colonies_' num2str(i)]);
DAPI = colonies(coli).loadImage(colDir, DAPIChannel);
imshow(DAPI,[]);

%% visualize moments for excluding asymmetric ones

ci = 4;
for coli = 15%:numel(colonies)

    colDir = fullfile(dataDir,['colonies_' num2str(i)]);
    img = colonies(coli).loadImage(colDir);
    colonies(coli).calculateMoments(img);

    CM = colonies(coli).CM{ci}; 
    I = colonies(coli).I{ci};
    [V,~] = eig(I);
    
    imshow(img(:,:,ci),[])
    hold on
    C = round([size(img,1) size(img,2)]/2);
    quiver(C(1),C(2),CM(1),CM(2),1,'g','LineWidth',2)
    quiver(C(1) + CM(1),C(1) + CM(2),V(1,1),V(2,1),100,'r','LineWidth',2)
    quiver(C(1) + CM(1),C(1) + CM(2),V(1,2),V(2,2),100,'r','LineWidth',2)
    hold off
end    
