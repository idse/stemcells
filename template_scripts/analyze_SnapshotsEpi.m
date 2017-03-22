clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

dataDir = '/Users/Idse/data_tmp/cellfate/pulses/170222_SyringePumpAendo2';
filelist = {'Process_1918.vsi','Process_1919.vsi','Process_1920.vsi','Process_1921.vsi'};

vsifile = fullfile(dataDir,filelist{1});
[~,barefname,~] = fileparts(vsifile);

% TODO : local normalization by DAPI

% metadata
%----------------

meta = Metadata(vsifile);
    
% manually entered metadata
%------------------------------

meta.channelLabel = {'Sox17','H2B','Smad4','Sox2'};
meta.conditions = {'mtesr','+A10', '+A100', '+A10-1hr-pulses'};

nucChannel = 2;
dataChannels = [2 1 4];

markerChannels = setdiff(dataChannels, nucChannel);

% set manual intensity limits so images can be compared
Ilim = {[100 2^16-1],[100 1500],[100 2000],[100 500]};

%% load data for parameter check

n = 1;
m = 2;
ymin = n*2048;
xmin = n*2048;
ymax = ymin + m*2048;
xmax = xmin + m*2048;
preview = {};
preview8bit = {};

series = 1;
img_bf = bfopen_mod(vsifile,xmin,ymin,xmax-xmin+1,ymax-ymin+1,series);

for ci = 1:numel(dataChannels)
    
    preview{ci} = img_bf{1}{dataChannels(ci),1};                    
    preview8bit{ci} = uint8((2^8-1)*mat2gray(preview{ci}, Ilim{dataChannels(ci)}));
end
filename = fullfile(dataDir, 'preview', [barefname '_RGBpreview.tif']);
imwrite(cat(3,preview8bit{:}), filename);

imshow(cat(3,preview8bit{:}))

%% load segmentation to test parameters

fi = 1;
vsifile = fullfile(dataDir, filelist{fi});
P = Position(meta.nChannels, vsifile, meta.nTime);
P.setID(pi);
seg = P.loadSegmentation(fullfile(dataDir,'MIP'), nucChannel);
segpart = seg(xmin:xmax, ymin:ymax);

%% finding the right options for extract data

% externally I will have all indices starting at 1
% the Andor offset to start at 0 will be internal

opts = struct(  'cytoplasmicLevels',    true,... 
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
figure, imshow(cat(3,   nucim + s*mat2gray(segpartclean),...
                        nucim,...
                        nucim + s*segpart(:,:,time)))

%% process data for first image

debugInfo = P.extractData(dataDir, nucChannel, opts);

nucmaskpart = debugInfo.nucmask(xmin:xmax, ymin:ymax);
figure, imshow(cat(3,   nucim + s*nucmaskpart,...
                        nucim,...
                        nucim + s*segpart(:,:,time)))
                    
% save
save(fullfile(dataDir, [barefname '_positions.mat']),'P');

%% scatter cell locations on top of part of preview to check result

XY = P.cellData.XY;

inpreview =     XY(:,1) < xmax & XY(:,1) > xmin...
                & XY(:,2) < ymax & XY(:,2) > ymin;
XY = XY(inpreview,:);

imshow(nucim,[])
hold on
scatter(XY(:,1) - xmin,XY(:,2) - ymin, '.r');
hold off

%% process the other images

for fi = 1:numel(filelist)
    
    vsifile = fullfile(dataDir, filelist{fi});
    [~,barefname,~] = fileparts(vsifile);

    P = Position(meta.nChannels, vsifile, meta.nTime);
    P.setID(pi);
    P.extractData(dataDir, nucChannel, opts);

    save(fullfile(dataDir, [barefname '_positions.mat']),'P');
end

%% load processed data for all files in the file list

allData = {};
nucLevelAll = [];

for i = 1:numel(filelist)

    vsifile = fullfile(dataDir,filelist{i});
    [~,barefname,~] = fileparts(vsifile);

    load(fullfile(dataDir, [barefname '_positions.mat']));
    allData{i} = P;
    
    N = P.cellData.nucLevel(:,nucChannel);
    nucLevel = bsxfun(@rdivide, P.cellData.nucLevel, N);
    nucLevelAll = cat(1, nucLevelAll, nucLevel);
end

% determine limits for displaying intensities
tol = 0.001;
lim = {};
for i = 1:meta.nChannels
    A = nucLevelAll(:,i);
    lim{i} = stretchlim(mat2gray(A),tol)*(max(A) - min(A)) + min(A);
end

%% histogram of different channels
% for w2_E6-CHIR-Act.vsi we see low nanog levels but many saturating sox17 cells

for fi = 1:numel(filelist)
    
    vsifile = fullfile(dataDir,filelist{fi});
    [~,barefname,~] = fileparts(vsifile);
    
    p = allData{fi};    
    N = p.cellData.nucLevel(:,nucChannel);
    nucLevel = bsxfun(@rdivide, p.cellData.nucLevel, N);

    for channelIndex = markerChannels

        bins = linspace(0,lim{channelIndex}(2),40);
        n = histc(nucLevel(:,channelIndex), bins);
        bar(bins,n)
        axis([bins(1) bins(end) 0 max(n)]);

    %    hist(p.cellData.nucLevel(:,channelIndex), 40)
        title(meta.channelLabel{channelIndex})
        filename = fullfile(dataDir, [barefname '_distribution_' meta.channelLabel{channelIndex}]);
        saveas(gcf, filename);
        saveas(gcf, [filename '.png']);
    end
end

%% scatter plots of markers against each other

%markerChannels = setdiff(dataChannels, nucChannel);

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
