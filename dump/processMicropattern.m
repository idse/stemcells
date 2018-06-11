clear all; close all;

%addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

%dataDir = 'img798';
dataDir = '/Volumes/IdseData/160318_micropattern_lefty/control';

colDir = fullfile(dataDir,'colonies');

btfname = fullfile(dataDir, 'Col1n2-797.btf');
%btfname = fullfile(dataDir, '160318_control.btf');

vsifile = fullfile(dataDir,'Image797.vsi');

metaDataFile = fullfile(dataDir,'metaData.mat');

% size limit of data chunks read in
maxMemoryGB = 4;

% TODO : local normalization by DAPI
% TODO : background subtraction by large scale Fourier transform or
% TODO : make objects save as struct to make sure we don't lose data       

%%
% metadata
%----------------

% read metadata
if exist(metaDataFile,'file')
    disp('loading previously stored metadata');
    load(metaDataFile);
elseif exist('vsifile','var') && exist(vsifile,'file')
    disp('extracting metadata from vsi file');
    meta = MetadataMicropattern(vsifile);
    
% or pretty much enter it by hand
else
    meta = MetadataMicropattern();
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

meta.channelLabel = {'DAPI','Cdx2','Sox2','Bra'};

meta.colRadiiMicron = 700;
meta.colMargin = 10; % margin outside colony to process, in pixels
meta.colRadiiPixel = meta.colRadiiMicron/meta.xres;

DAPIChannel = find(strcmp(meta.channelNames,'DAPI'));

save(fullfile(dataDir,'metaData'),'meta');

%% processing loop

% masks for radial averages
[radialMaskStack, edges] = makeRadialBinningMasks(meta);

% split the image up in big chunks for efficiency
maxBytes = (maxMemoryGB*1024^3);
bytesPerPixel = 2;
dataSize = meta.ySize*meta.xSize*meta.nChannels*bytesPerPixel;
nChunks = ceil(dataSize/maxBytes);

if nChunks > 1
    nRows = 2;
else
    nRows = 1;
end
nCols = ceil(nChunks/nRows);

xedge = (0:nCols)*(meta.xSize/nCols);
yedge = (0:nRows)*(meta.ySize/nRows);

% define the data structures to be filled in
preview = zeros(floor([2048 2048*meta.xSize/meta.ySize 4]));
mask = false([meta.ySize meta.xSize]);
colonies = [];
chunkColonies = {};
% L = 2048;
% bg = (2^16-1)*ones([L L meta.nChannels],'uint16');

chunkIdx = 0;

for n = 1:numel(yedge)-1
    for m = 1:numel(xedge)-1

        chunkIdx = chunkIdx + 1;

        disp(['reading chunk ' num2str(chunkIdx) ' of ' num2str(nChunks)])
        
        xmin = uint32(xedge(m) + 1); xmax = uint32(xedge(m+1));
        ymin = uint32(yedge(n) + 1); ymax = uint32(yedge(n+1));
        
        % add one sided overlap to make sure all colonies are completely
        % within at least one chunk
        % theoretically one radius should be enough but we need a little
        % margin
        if n < nRows 
            ymax = ymax + 1.25*max(meta.colRadiiPixel); 
        end
        if m < nCols
            xmax = xmax + 1.25*max(meta.colRadiiPixel);
        end
        chunkheight = ymax - ymin + 1;
        chunkwidth = xmax - xmin + 1;
        
        % for preview (thumbnail)
        ymaxprev = ceil(size(preview,1)*double(ymax)/meta.ySize);
        yminprev = ceil(size(preview,1)*double(ymin)/meta.ySize);
        xmaxprev = ceil(size(preview,2)*double(xmax)/meta.xSize);
        xminprev = ceil(size(preview,2)*double(xmin)/meta.xSize);
        
        img = zeros([chunkheight, chunkwidth, meta.nChannels],'uint16');
        for ci = 1:meta.nChannels
            tic
            disp(['reading channel ' num2str(ci)])
            img(:,:,ci) = imread(btfname,'Index',ci,'PixelRegion',{[ymin,ymax],[xmin, xmax]});
            preview(yminprev:ymaxprev,xminprev:xmaxprev, ci) = ...
                imresize(img(:,:,ci),[ymaxprev-yminprev+1, xmaxprev-xminprev+1]);
            toc
        end
        
        % determine background
        % bg = getBackground(bg, img, L);

        disp('determine threshold');
        if m == 1 && n == 1
            
            for ci = DAPIChannel
                
                imsize = 2048;
                si = size(img);
                xx = min(si(1),4*imsize); yy = min(si(2),4*imsize);
                if xx < 2*imsize
                    x0 = 1;
                else
                    x0 = imsize + 1;
                end
                if yy < 2*imsize
                    y0 = 1;
                else
                    y0 = imsize + 1;
                end
                    
                forIlim = img(x0:xx,y0:yy,ci);
                minI = double(min(forIlim(:)));
                maxI = double(max(forIlim(:)));
                forIlim = mat2gray(forIlim);
                
                IlimDAPI = (stretchlim(forIlim)*maxI + minI)';
                t = 0.8*graythresh(forIlim)*maxI + minI;
            end
        end
        mask(ymin:ymax,xmin:xmax) = img(:,:,DAPIChannel) > 0.8*t;

        disp('find colonies');
        tic
        s = round(20/meta.xres);
        range = [xmin, xmax, ymin, ymax];
        [chunkColonies{chunkIdx}, cleanmask] = findColonies(mask, range, meta, s);
        toc

        disp('merge colonies')
        prevNColonies = numel(colonies);
        if prevNColonies > 0
            D = distmat(cat(1,colonies.center), chunkColonies{chunkIdx}.center);
            [i,j] = find(D < max(meta.colRadiiPixel)*meta.xres);
            chunkColonies{chunkIdx}(j) = [];
        end
        % add fields to enable concatenating
        colonies = cat(2,colonies,chunkColonies{chunkIdx});
        
        disp('process individual colonies')
        tic
        
        % channels to save to individual images
        if ~exist(colDir,'dir')
            mkdir(colDir);
        end

        nColonies = numel(colonies);
        
        for coli = prevNColonies+1:nColonies
            
            % store the ID so the colony object knows its position in the
            % array (used to then load the image etc)
            colonies(coli).setID(coli);
            
            fprintf('.');
            if mod(coli,60)==0
                fprintf('\n');
            end
            
            colDiamMicron = 2*colonies(coli).radiusMicron; 
            
            b = colonies(coli).boundingBox;
            colnucmask = mask(b(3):b(4),b(1):b(2));
            colcytmask = imdilate(colnucmask,strel('disk',round(5/meta.xres)))-colnucmask;
            
            b(1:2) = b(1:2) - double(xmin - 1);
            b(3:4) = b(3:4) - double(ymin - 1);
            colimg = img(b(3):b(4),b(1):b(2), :);
            
            colmaskClean = cleanmask(b(3):b(4),b(1):b(2));

            % write colony image
            colonies(coli).saveImage(colimg, colDir);

            % write DAPI separately for Ilastik
            colonies(coli).saveImage(colimg, colDir, DAPIChannel);

            % do radial binning
            colType = find(meta.colRadiiMicron == colonies(coli).radiusMicron);
            N = size(radialMaskStack{colType},3);
            nucradavg = zeros([N meta.nChannels]);
            nucradstd = zeros([N meta.nChannels]);
            cytradavg = zeros([N meta.nChannels]);
            cytradstd = zeros([N meta.nChannels]);
            
            % store bin edges, to be reused by segmented profiles later
            colonies(coli).radialProfile.BinEdges = edges{colType};
            
            for ri = 1:N
                % for some reason linear indexing is faster than binary
                colnucbinmask = find(radialMaskStack{colType}(:,:,ri) & colnucmask);
                colcytbinmask = find(radialMaskStack{colType}(:,:,ri) & colcytmask);
                
                for ci = 1:meta.nChannels
                    imc = colimg(:,:,ci);
                    % most primitive background subtraction: minimal value
                    % within the colony
                    % min(imc(colmaskClean)) doubles the computatation time
                    imc = imc - min(imc(:));
                    
                    imcbin = imc(colnucbinmask);
                    nucradavg(ri,ci) = mean(imcbin);
                    nucradstd(ri,ci) = std(double(imcbin));
                    
                    imcbin = imc(colcytbinmask);
                    cytradavg(ri,ci) = mean(imcbin);
                    cytradstd(ri,ci) = std(double(imcbin));
                end
            end
            colonies(coli).radialProfile.NucAvg = nucradavg;
            colonies(coli).radialProfile.NucStd = nucradstd;
            colonies(coli).radialProfile.CytAvg = cytradavg;
            colonies(coli).radialProfile.CytStd = cytradstd;
        end
        fprintf('\n');
        toc
    end
end

save(fullfile(dataDir,'colonies'), 'colonies');

%% visualize
    
figure,
imshow(imadjust(mat2gray(preview(:,:,DAPIChannel))))
hold on
CM = cat(1,colonies.center);
scale = size(preview,1)/meta.ySize;
CM(:,2) = CM(:,2)*scale;
CM(:,1) = CM(:,1)*scale;
radius = cat(1,colonies.radiusPixel)*scale;    
%scatter(CM(:,1),CM(:,2),'.r');
viscircles(CM,radius,'LineWidth',1)
hold off
saveas(gcf,fullfile(dataDir,'coloniesOverviewDAPI.png'));
hold on 
for i = 1:size(CM,1)
   text(CM(i,1),CM(i,2),num2str(i),'Color','red','BackgroundColor','white',...
       'Margin',1,'FontSize',5,'HorizontalAlignment','center'); 
end
hold off
saveas(gcf,fullfile(dataDir,'coloniesOverviewDAPIlabels.png'));

%% visualize
    
figure,
previewRGB = preview(:,:,2:4);
for i = 1:3
    previewRGB(:,:,i) = imadjust(mat2gray(previewRGB(:,:,i)));
end
channelPermutation = [3 1 2];
imshow(previewRGB(:,:,channelPermutation))
hold on
CM = cat(1,colonies.center);
scale = size(preview,1)/meta.ySize;
CM(:,2) = CM(:,2)*scale;
CM(:,1) = CM(:,1)*scale;
hold off
d = 80;
text(1,1,meta.channelLabel{channelPermutation(1)+1},'Color','r','VerticalAlignment','top','BackgroundColor','white')
text(1,d,meta.channelLabel{channelPermutation(2)+1},'Color','g','VerticalAlignment','top','BackgroundColor','white')
text(1,2*d,meta.channelLabel{channelPermutation(3)+1},'Color','b','VerticalAlignment','top','BackgroundColor','white')
saveas(gcf,fullfile(dataDir,'coloniesOverviewRGB.png'));

%% also save preview composite

preview = uint16(preview);
fname = fullfile(dataDir,'coloniesOverview.tif');
imwrite(preview(:,:,1),fname);
for i = 2:4
    imwrite(preview(:,:,i),fname,'WriteMode','Append');
end

%% read segmentation
% 
% make data in classifier 'copied to protect file'

if ~exist('colonies','var')
    load(fullfile(dataDir,'colonies'));
end

%%
% take just the large colonies

colRadii = cat(1,colonies.radiusMicron);
colonies1000idx = colRadii == 500;
colonies1000 = colonies(colonies1000idx);

%%
coli = 4;

DAPI = colonies(coli).loadImage(colDir, DAPIChannel);
seg = colonies(coli).loadSegmentation(colDir, DAPIChannel);

figure, imshow(label2rgb(bwlabel(seg),'jet','k','shuffle'));

%%
options = struct('minAreaStd', 1, 'minSolidity',0);
seg2 = separateFusedNuclei(seg, options);
figure, imshow(label2rgb(bwlabel(seg2),'jet','k','shuffle'));

%% extract data using Ilastik segmentation

tic

cleanupOpts = struct('minArea', 30, 'cytoplasmicLevels', false);
opts = struct('cleanupOptions', cleanupOpts);

for coli = [colonies1000.ID]
    
    % extract segmented data
    colonies(coli).extractData(colDir, DAPIChannel, opts);

    % radial binning of segmented data
    colType = find(meta.colRadiiMicron == colonies(coli).radiusMicron);
    colonies(coli).makeRadialAvgSeg()
end

save(fullfile(dataDir,'colonies'), 'colonies');

%% display segmented radial profile

colRadii = cat(1,colonies.radiusMicron);
colonies1000idx = colRadii == 500;
colonies1000 = colonies(colonies1000idx);

i = find(meta.colRadiiMicron == 500);
r = imfilter(colonies1000(1).radialProfile.BinEdges,[1 1]/2)*meta.xres;
r(1) = colonies1000(1).radialProfile.BinEdges(1)*meta.xres;
r = r(1:end-1);

colCat = cat(3,colonies1000(:).radialProfile);
colCat = cat(3,colCat.NucAvgSeg);
avg = mean(colCat,3);
avgNormalizedSeg = bsxfun(@rdivide, avg, avg(:,1));

plot(r, avgNormalizedSeg(:,2:4))
legend(meta.channelLabel(2:4));
axis([min(r) max(r) 0 4]);

%% not segmented

figure,
colCat = cat(3,colonies1000(:).radialProfile);

nucAvgAll = mean(cat(3,colCat.NucAvg),3);
cytAvgAll = mean(cat(3,colCat.CytAvg),3);

nucAvgAllNormalized = bsxfun(@rdivide, nucAvgAll, nucAvgAll(:,DAPIChannel));

% how normalize cytoplasmic levels?
cytAvgAllNormalized = bsxfun(@rdivide, cytAvgAll, nucAvgAll(:,DAPIChannel));

plot(r, nucAvgAllNormalized(:,2:4))
legend(meta.channelLabel(2:4));
axis([min(r) max(r) 0 2]);
hold on
plot(r, cytAvgAllNormalized(:,2:4),'--')
hold off

%% compare

figure,
plot(r, avgNormalizedSeg(:,2:4))
axis([min(r) max(r) 0 2]);
hold on
plot(r, nucAvgAllNormalized(:,2:4),'--')
hold off
legend(meta.channelLabel(2:4));
title('segmented vs not segmented (dashed)')

%% look at single segmented colony

coli = 262;
avgSeg = colonies(coli).radialProfile.NucAvgSeg;
avgSegNormalized = bsxfun(@rdivide, avgSeg, avgSeg(:,1));
figure(3)
plot(r, avgSegNormalized(:,2:4))
legend(meta.channelLabel(2:4))
axis([min(r) max(r) 0 3]);

%% compare non-segmented

avg = colonies(coli).radialProfile.NucAvg;
nucAvgAllNormalized = bsxfun(@rdivide, avg, avg(:,1));
figure(2)
plot(r, nucAvgAllNormalized(:,2:4))
legend(meta.channelLabel(2:4))
axis([min(r) max(r) 0 3]);
