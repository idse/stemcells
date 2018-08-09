clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

dataDir = '/Users/idse/data_tmp/0_George/180702_IWP2 Quants';
filenrs = [4716 4719 4722 4725 4728 4731];

conditions = strcat({'72','66','60','48','36','24'},'h IWP2');

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
doubleNormalize = true;

for i = 1:numel(filenrs)
    load(fullfile(dataDir,['colonies_' num2str(filenrs(i))]), 'colonies');

    % make mean radial profiles
    figure,
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

load(fullfile(dataDir,'metaData.mat'),'meta');
doubleNormalize = true;

for i = 1:numel(filenrs)    
    load(fullfile(dataDir,['colonies_' num2str(filenrs(i))]), 'colonies');
    coloniesCombined{i} = colonies;
end

figure('Position',[0 0 1600 300]),
plotMultipleAveragesNoSegmentation(meta, colRadius, DAPIChannel,...
                                coloniesCombined, conditions, doubleNormalize)
saveas(gcf,fullfile(dataDir,'radialProfilesCombined.png'));

%% load colonies and show side by side

load(fullfile(dataDir,'metaData.mat'),'meta');

% representative colonies for each condition
repcols = [7 4 6 6 7 7; 1 1 1 1 1 1; 12 12 12 12 12 12]'; 
ci = 1;

figure('Position',[0 0 1600 1000]),
m = 3;
n = numel(filenrs);

for j = 1:m
    for i = 1:n

        ii = (j-1)*n + i;

        load(fullfile(dataDir,['colonies_' num2str(filenrs(i)) '.mat'])); 
        colDir = fullfile(dataDir,['colonies_' num2str(filenrs(i))]);

        coli = repcols(ii);
        img = colonies(coli).loadImage(colDir, ci);
        bg = imopen(img, strel('disk',15));
        img = img - bg;
        if j==1 && i == 1
            Ilim = round(stretchlim(img)*(2^16-1));
        end
        subplot_tight(m,n,ii)
        imshow(img,Ilim);
        title(conditions{i})
    end
end
saveas(gcf,fullfile(dataDir,['compareColoniesImopenbg_' meta.channelLabel{ci} '.png']));

%% try background subtraction

i = 1;
load(fullfile(dataDir,['colonies_' num2str(filenrs(i)) '.mat'])); 

coli = 7;
colDir = fullfile(dataDir,['colonies_' num2str(filenrs(i))]);
figure,
im = colonies(coli).loadImage(colDir, 3);
im = mat2gray(im);

bg = imopen(im,strel('disk',15));
imshow(im-bg,[])

% sigma = 10;
% bg = imfilter(im,fspecial('gauss',3*sigma, sigma));
% bgsub = im-bg;
% %bgsub(bgsub < 0) = 0;
% imshow(bgsub,[])

% sigma = 3;
% bg = imfilter(im,fspecial('log',3*sigma, sigma));
% imshow(imadjust(mat2gray(im - 100*bg)),[]);

%imshow(imadjust(mat2gray(pixelNorm)),[])


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
