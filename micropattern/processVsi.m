function processVsi(vsifile,dataDir,varargin)
% Function that takes a vsifile as input and outputs to dataDir.
% will create a colonies.mat file with the colonies array, a metaData.mat
% file with variable meta and a directory colonies containing images of
% individual colonies.
% varargin allows for overwrite of default metadata options. Must have keyword:value pairs.
%   -channelLabel (default {'DAPI','Cdx2','Sox2','Bra'})
%   -colRadiiMicron (default [200 500 800 1000]/2 )
%   -colMargin (default 10)
%   -DAPI channel - nuclear marker channel. By default looks for DAPI in
%   channel labels
%   -metaDataFile: can read meta structure from file and skip extracting
%   from vsi. (Default reads from vsi).

%%
% metadata
%----------------

% read metadata inputs
in_struct = varargin2parameter(varargin);

if isfield(in_struct,'metaDataFile')
    disp('loading previously stored metadata');
    metaDataFile = in_struct.metaDataFile;
    load(metaDataFile);
    metadatadone = true;
    meta2 = MetadataMicropattern(vsifile);
    %make sure image size is correct if importing metadata
    meta.xSize = meta2.xSize;
    meta.ySize = meta2.ySize;
else
    disp('extracting metadata from vsi file');
    meta = MetadataMicropattern(vsifile);
    metaDataFile = fullfile(dataDir,'metaData.mat');
    metadatadone = false;
end
% or pretty much enter it by hand

% manually entered metadata
%------------------------------

if ~metadatadone
    %defaults
    meta.channelLabel = {'DAPI','Cdx2','Sox2','Bra'};
    meta.colRadiiMicron = [200 500 800 1000]/2;
    meta.colMargin = 10; % margin outside colony to process, in pixels
    
    %override from inputs
    if isfield(in_struct,'channelLabel')
        meta.channelLabel = in_struct.channelLabel;
    end
    if isfield(in_struct,'colRadiiMicron')
        meta.colRadiiMicron = in_struct.colRadiiMicron;
    end
    if isfield(in_struct,'colMargin')
        meta.colMargin = in_struct.colMargin;
    end
    
    if isfield(in_struct,'channelNames')
        meta.channelNames = in_struct.channelNames;
        meta.nChannels = length(meta.channelNames);
    end
    
    meta.colRadiiPixel = meta.colRadiiMicron/meta.xres;
    
end
if ~exist(dataDir)
    mkdir(dataDir);
end
save(fullfile(dataDir,'metaData.mat'),'meta');

maxMemoryGB = 4;

if isfield(in_struct,'DAPIChannel')
    DAPIChannel = in_struct.DAPIChannel;
else
    DAPIChannel = find(strcmp(meta.channelNames,'DAPI'));
end
colDir = fullfile(dataDir,'colonies');
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
prevNWells = 0;
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
        [~, ext] = strtok(vsifile,'.');
        if strcmp(ext,'.vsi')
            img_bf = bfopen_mod(vsifile,xmin,ymin,xmax-xmin+1,ymax-ymin+1,1); %read only the 1st series from the vsi
            
            for ci = 1:meta.nChannels
                %             tic
                %             disp(['reading channel ' num2str(ci)])
                %             img(:,:,ci) = imread(btfname,'Index',ci,'PixelRegion',{[ymin,ymax],[xmin, xmax]});
                img(:,:,ci) = img_bf{1}{ci,1};
                preview(yminprev:ymaxprev,xminprev:xmaxprev, ci) = ...
                    imresize(img(:,:,ci),[ymaxprev-yminprev+1, xmaxprev-xminprev+1]);
                %            toc
            end
        else
            for ci = 1:meta.nChannels
                img(:,:,ci) = imread(vsifile,ci);
                preview(yminprev:ymaxprev,xminprev:xmaxprev, ci) = ...
                    imresize(img(:,:,ci),[ymaxprev-yminprev+1, xmaxprev-xminprev+1]);
            end
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
        [chunkColonies{chunkIdx}, cleanmask, welllabel] = findColonies(mask, range, meta, s);
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
            
            
            %increment well to account for previous
            colonies(coli).well = colonies(coli).well + prevNWells;
            
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
        prevNWells = prevNWells+max(max(welllabel));
        
        fprintf('\n');
        toc
        
    end
end
preview = uint16(preview);
imwrite(squeeze(preview(:,:,1)),fullfile(dataDir,'previewDAPI.tif'));
imwrite(preview(:,:,2:4),fullfile(dataDir,'previewRGB.tif'));

save(fullfile(dataDir,'colonies'), 'colonies');