clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

dataDir = '/Users/Idse/data_tmp/cellfate/pulses/170222_SyringePumpAendo2';
filelist = {'Process_1918.vsi','Process_1919.vsi','Process_1920.vsi','Process_1921.vsi'};

vsifile = fullfile(dataDir,filelist{1});
[~,barefname,~] = fileparts(vsifile);

resultsDir = fullfile(dataDir,'results');
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end
if ~exist(fullfile(resultsDir,'figs'),'dir')
    mkdir(fullfile(resultsDir,'figs'));
end

% TODO : local normalization by DAPI

% metadata
%----------------

meta = Metadata(vsifile);
    
% manually entered metadata
%------------------------------

meta.channelLabel = {'Sox17','H2B','Smad4','Sox2'};
meta.conditions = {'mtesr','+A10', '+A100', '+A10-1hr-pulses'};

nucChannel = 2;
dataChannels = [1 4 2];

markerChannels = setdiff(dataChannels, nucChannel);

%% load data for parameter check

n = 1;
m = 2;
ymin = n*2048;
xmin = n*2048;
ymax = ymin + m*2048;
xmax = xmin + m*2048;
preview = {};
Ilim = {};

series = 1;

% load data and determine intensity limits per file
for fi = 1:numel(meta.conditions)

    vsifile = fullfile(dataDir,filelist{fi});
    img_bf = bfopen_mod(vsifile,xmin,ymin,xmax-xmin+1,ymax-ymin+1,series);
    
    for ci = 1:numel(dataChannels)

        im = img_bf{1}{dataChannels(ci),1};
        preview{ci, fi} = im;
        maxim = double(max(im(:)));
        minim = double(min(im(:)));
        Ilim{dataChannels(ci), fi} = stretchlim(mat2gray(im))*(maxim-minim) + minim;
    end
end

%% combine lookup tables and save RGB previews

preview8bit = {}; 
for fi = 1:numel(meta.conditions)
    
    for ci = 1:numel(dataChannels)  
        
        Imin = min(min([Ilim{dataChannels(ci),:}]));
        Imax = max(max([Ilim{dataChannels(ci),:}]));
        preview8bit{ci,fi} = uint8((2^8-1)*mat2gray(preview{ci, fi}, [Imin Imax]));
    end

    previewChannels = 1:3;
    [~,barefname,~] = fileparts(filelist{fi});
    filename = fullfile(dataDir, 'preview', ['RGBpreview_' barefname '_' [meta.channelLabel{dataChannels}] '.tif']);
    RGBim = cat(3,preview8bit{previewChannels,fi});
    imwrite(RGBim, filename);
end

%% smaller previews to copy to slide

n = 1; m = numel(meta.conditions);
L = 2^12;
screensize = get( 0, 'Screensize' );
margin = 50;
fs = 15;
figure('Position', [1, 1, screensize(3), screensize(3)/m + margin/2]);
for fi = 1:numel(meta.conditions)
    
    subplot_tight(n,m,fi)
    RGBim = cat(3,preview8bit{previewChannels,fi});
    imshow(RGBim(1:L,1:L,:));
    title(meta.conditions(fi),'FontSize',fs,'FontWeight','bold');
    labelstr = ['\color{red}'   meta.channelLabel{dataChannels(1)}...
                '\color{green}' meta.channelLabel{dataChannels(2)}...
                '\color[rgb]{0.1, 0.5, 1}' meta.channelLabel{dataChannels(3)}];
    text(margin, L - 2*margin, labelstr,'FontSize',fs,'FontWeight','bold');
end
saveas(gcf, fullfile(resultsDir,'overview.png'));
close;

%% load segmentation to test parameters

vsifile = fullfile(dataDir, filelist{fi});
P = Position(meta.nChannels, vsifile, meta.nTime);
P.setID(pi);
seg = P.loadSegmentation(fullfile(dataDir,'MIP'), nucChannel);
segpart = seg(xmin:xmax, ymin:ymax);

%% finding the right options for extract data

% externally I will have all indices starting at 1
% the Andor offset to start at 0 will be internal

opts = struct(  'cytoplasmicLevels',    false,... 
                    'dataChannels',     1:meta.nChannels,... 
                    'tMax',             meta.nTime,...
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'nucShrinkage',     2,...
                    'cytoSize',         5,...
                    'cytoMargin',       5,...
                    'bgMargin',         4);
                
opts.cleanupOptions = struct('separateFused', true,'openSize',5,...
    'clearBorder',true, 'minAreaStd', 1, 'minSolidity',0, 'minArea',500);

% try out the nuclear cleanup settings on some frame:
segpartclean = nuclearCleanup(segpart, opts.cleanupOptions);

% visualize
s = 0.2;
nucim = mat2gray(preview8bit{dataChannels==nucChannel});

RGBsegCheck = cat(3,   nucim + s*mat2gray(segpartclean),...
                        nucim,...
                        nucim + s*segpart);
figure, imshow(RGBsegCheck)

filename = fullfile(resultsDir, [barefname 'segCheck.tif']);
imwrite(RGBsegCheck, filename);

%% process data for first image

debugInfo = P.extractData(dataDir, nucChannel, opts);

nucmaskpart = debugInfo.nucmask(xmin:xmax, ymin:ymax);
figure, imshow(cat(3,   nucim + s*nucmaskpart,...
                        nucim,...
                        nucim + s*segpart))
                    
% save
save(fullfile(resultsDir, [barefname '_positions.mat']),'P');

%% scatter cell locations on top of part of preview to check result

XY = P.cellData.XY;

inpreview =     XY(:,1) < xmax & XY(:,1) > xmin...
                & XY(:,2) < ymax & XY(:,2) > ymin;
XY = XY(inpreview,:);

imshow(nucim,[])
hold on
scatter(XY(:,1) - xmin,XY(:,2) - ymin, 200, '.r');
hold off

%% process the other images

for fi = 2:numel(filelist)
    
    vsifile = fullfile(dataDir, filelist{fi});
    [~,barefname,~] = fileparts(vsifile);

    P = Position(meta.nChannels, vsifile, meta.nTime);
    P.setID(pi);
    P.extractData(dataDir, nucChannel, opts);

    save(fullfile(resultsDir, [barefname '_positions.mat']),'P');
end

%% make distributions

% load all data  
allData = {};
for fi = 1:numel(filelist)

    vsifile = fullfile(dataDir, filelist{fi});
    [~,barefname,~] = fileparts(vsifile);
    load(fullfile(dataDir,'results', [barefname '_positions.mat']));
    allData{fi} = P;
end

% make distributions out of processed data
stats = cellStats(allData, meta, markerChannels);
stats.makeHistograms(40);

% plot distributions

% overlay of the distribution of different conditions
for channelIndex = markerChannels

    stats.plotDistributionComparison(channelIndex)

    filename = ['distOverlay_' meta.channelLabel{channelIndex}];
    saveas(gcf, fullfile(resultsDir, 'figs', filename));
    saveas(gcf, fullfile(resultsDir, [filename '.png']));
    close;
    
    cumulative = true;
    stats.plotDistributionComparison(channelIndex, cumulative)

    filename = ['cumdistOverlay_' meta.channelLabel{channelIndex}];
    saveas(gcf, fullfile(resultsDir, 'figs', filename));
    saveas(gcf, fullfile(resultsDir, [filename '.png']));
    close;
end

%% try clustering in 2D

k = 2;
regularizationValue = 10^(-3);
cutoff = 0.95;
stats.makeClusters(k, regularizationValue, cutoff);

figure, 

for ci1 = 1:numel(markerChannels)
    for ci2 = ci1+1:numel(markerChannels)

        clf 
        c1 = markerChannels(ci1);
        c2 = markerChannels(ci2);
        showClusters = true;
        for conditionIdx = 1:numel(stats.conditions)
            stats.makeScatterPlot(conditionIdx, [c1 c2], showClusters)
        end
        title('all data clustered');

        [~,barefname,~] = fileparts(filelist{conditionIdx});
        filename = ['cluster_' meta.channelLabel{c1} '_' meta.channelLabel{c2}];
        saveas(gcf, fullfile(resultsDir, 'figs', filename));
        saveas(gcf, fullfile(resultsDir, [filename '.png']));
        close;
    end
end

%% if good visualize clusters for each condition

figure,

for ci1 = 1:numel(markerChannels)
    for ci2 = ci1+1:numel(markerChannels)
        
        channelIdx = markerChannels([ci1 ci2]);
        showClusters = true;
        
        for conditionIdx = 1:numel(stats.conditions)

            clf
            stats.makeScatterPlot(conditionIdx, channelIdx, showClusters)
            
            [~,barefname,~] = fileparts(filelist{conditionIdx});
            filename = ['cluster_' stats.channelLabel{c1} '_' stats.channelLabel{c2} '_' barefname];
            saveas(gcf, fullfile(resultsDir, 'figs', filename));
            saveas(gcf, fullfile(resultsDir, [filename '.png']));
            close;
        end
    end
end

% plot distributions per cluster
stats.makeClusterHistograms();

for clusterLabel = 1:stats.nClusters
    for channelIndex = markerChannels

        cumulative = false;
        stats.plotDistributionComparison(channelIndex, cumulative, clusterLabel)

        filename = ['distOverlay_' meta.channelLabel{channelIndex} '_cluster' num2str(clusterLabel)];
        saveas(gcf, fullfile(resultsDir, 'figs', filename));
        saveas(gcf, fullfile(resultsDir, [filename '.png']));
        close;

        cumulative = true;
        stats.plotDistributionComparison(channelIndex, cumulative, clusterLabel)

        filename = ['cumdistOverlay_' meta.channelLabel{channelIndex} '_cluster' num2str(clusterLabel)];
        saveas(gcf, fullfile(resultsDir, 'figs', filename));
        saveas(gcf, fullfile(resultsDir, [filename '.png']));
        close;
    end
end

% summary 
txtfile = fullfile(resultsDir, 'statsSummary.txt');
if exist(txtfile,'file')
    delete(txtfile);
end
diary(txtfile);
diary on
stats.printSummary();
diary off
