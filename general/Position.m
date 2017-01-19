classdef Position < handle
    % Data class to store cell data in a field of view

    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------

    properties

        ncells              % number of cells
        density             % cell density

        cellData            % structure array indexed by time:
                            % - XY
                            % - area
                            % - nucLevel:       nCell x nChannels array
                            % - cytLevel
                            % - background:     nChannels long vector 
                            % - nucLevelAvg
                            % - cytLevelAvg
                            
        timeTraces          % structure: reorganizes cellData 
                            % - nucLevelAvg:    vector indexed by time
                            % - cytLevelAvg 
                            % - background 
                            % NOT IMPLEMENTED YET: traces of tracked cells
                            
        dataChannels        % we may not want to quantify all channels

        % different from cellTracker: 
        %------------------------------

        nChannels           % number of colors imaged
        nTime               % number of time points
        filename   
    end

    properties (SetAccess = protected)
        ID                  % identifier of colony 
    end

    properties (Dependent)
        
        bareFilename;       % filename without extension
        %data                % cell by cell data for colony
                            % data cols: x, y, area, -1, 
                            % (nuclear, cytoplasmic) mean for each channel
    end

    methods

        % constructor
        function this = Position(nChannels, filename, nTime)
            % Position(nChannels, filename, nTime)
            
            % matlab sucks
            if nargin == 0
                return
            end
            
            if ~exist('nTime','var')
                nTime = 1;
                this.timeTraces = [];
            else
                this.timeTraces = struct();
            end

            if ~isnumeric(nChannels)
                error('first argument is nChannels, which is a number');
            end
            if ~isnumeric(nTime)
                error('third argument is nTime, which is a number');
            end
            
            this.nChannels = nChannels;
            this.nTime = nTime;
            
            % strip off path in case it was provided
            [~, name, ext] = fileparts(filename);
            this.filename = [name ext];
            
            this.cellData = struct();
        end

        % saving and loading
        %---------------------------------

        function img = loadImage(this, dataDir, channels, time)
            % load position image
            %
            % img = loadImage(dataDir, channels)
            % img = loadImage(dataDir, channels, time)
            %
            % dataDir:  main data directory 
            % channels: desired channels to be loaded, leave empty for all
            % time:     time to load
            %
            % img:      loaded image

            if ~exist('channels','var')
                channels = 1:this.nChannels;
            end
            fname = fullfile(dataDir, this.filename);
            [~,~,ext] = fileparts(this.filename);
            
            if strcmp(ext,'.tif') || strcmp(ext,'.btf')

                info = imfinfo(fname);
                w = info.Width;
                h = info.Height;

                img = zeros([h w numel(channels)],'uint16');
                for cii = 1:numel(channels)
                    img(:,:,cii) = imread(fname, channels(cii));
                end
                
                if exist('time','var') && time > 1
                    error('todo : include reading for dynamic not Andor or epi');
                end
                
            elseif strcmp(ext,'.vsi')
                
                r = bfGetReader(fname);
                img = zeros([r.getSizeY() r.getSizeX() numel(channels) r.getSizeZ()], 'uint16');
                for cii = 1:numel(channels)
                    for zi = 1:r.getSizeZ()
                        img(:,:,cii,zi) = bfGetPlane(r, r.getIndex(zi-1,channels(cii),time-1));
                    end
                end
                r.close();
                img = squeeze(img);
% OLD, probably obsolete (reads vsi, but can't deal with time)
%             else
%                 %img = readStack(fname);
%                 xmin = []; ymin = []; xmax = []; ymax = [];
%                 series = 1;
%                 img_bf = bfopen_mod(fname,xmin,ymin,xmax-xmin+1,ymax-ymin+1,series,channels);
%                 img = cat(3,img_bf{1}{:,1});
            end
            
            %disp(['loaded image ' fname]);
        end

        function seg = loadSegmentation(this, dataDir, channel)
            % load segmentation
            %
            % seg = loadSegmentation(dataDir, channel)
            %
            % dataDir:  main data directory
            %
            % seg:      binary image of nuclei
            %
            % segmentation is expected to have same name as raw data
            % but with .tif -> '_Simple Segmentation.h5'
            %
            % by convention use green for foreground in Ilastik
            
            % for any normal filename this just takes of the extension
            % for a filename format i.e. bla_something_x%.4d_whatever.tif, 
            % it keeps bla_something

            if ~exist('channel','var')
                error('please provide channel');
            end
            
            [~,~,ext] = fileparts(this.filename);
            
            % just a convention, batchMIP_epi renames MIPs to Andor
            % convention with barefname the dataDir name, which is 
            % one level above /MIP
            if strcmp(ext,'.vsi')
                s = strsplit(dataDir,filesep);
                barefname = sprintf([s{end-1}, '_MIP_p%.4d'], this.ID-1); 
            else
                s = strsplit(this.filename,'_.%\.4d|\.','DelimiterType','RegularExpression');            
                barefname = s{1};
            end
            
            % for dynamic data the raw data will usually be a MIP
            % this code tries both
            listing = dir(fullfile(dataDir,[barefname '*h5']));
            if isempty(listing)
                s = strsplit(this.filename,'_.[0-9]+','DelimiterType','RegularExpression');
                barefname = sprintf([s{1} '_MIP_p%.4d'], this.ID-1);
                listing = dir(fullfile(dataDir,[barefname '_*h5']));
            end
            if isempty(listing)
                error(['segmentation ' [barefname '*h5'] ' for channel ' num2str(channel) ' not found in ' dataDir]);
            end

            seg = [];
            
            for i = 1:numel(listing)
                % DON'T COMMENT THE IF STATEMENT BELOW, IT WILL LOAD THE WRONG
                % CHANNEL
                if ~isempty(strfind(listing(i).name,sprintf('_w%.4d',channel-1)))... 
                        || ~isempty(strfind(listing(i).name,sprintf('_c%d',channel)))

                    fname = fullfile(dataDir,listing(i).name);
                    seg = h5read(fname, '/exported_data');
                    % Probabilities
                    if size(seg,1) > 1 
                        ssegsize = size(seg);
                        ssegsize(1) = 1;
                        simpleseg = false(ssegsize);
                        simpleseg(seg(2,:,:,:) >= 0.5) = true;
                        seg = squeeze(simpleseg);
                        
                    % Simple Segmentation
                    else                
                        seg = squeeze(seg == 2);
                    end
                end
            end
            
            if isempty(seg)
                warning(['segmentation for channel ' num2str(channel) ' not found in ' dataDir, ', may be naming convention problem']);
            else
                % Ilastik output data has xy transposed
                seg = permute(seg, [2 1 3]);
                disp(['loaded segmentation ' fname]);
            end
        end
        
        function MIPidx = loadMIPidx(this, dataDir, channel, time)
            % load MIP index of position
            %
            % MIPidx = loadMIPidx(this, dataDir, channel)
            %
            % dataDir:  main data directory 
            % channels: desired channels to be loaded, leave empty for all
            % time:     time to load, default 1 for static images
            %
            % img:      loaded image

            if ~exist('time','var')
                time = 1;
            end
            
            fname = this.filename;

            % remove time part of filename:
            [startIndex,endIndex]  = regexp(this.filename,'_t([0-9]+|%.4d)');
            if ~isempty(startIndex)
                fname = fname([1:startIndex-1 endIndex+1:end]);
            end
            
            % MIPidx is supposed to be single file for all times
            startIndex = regexp(fname,'_.[0-9]+');
            fname = [fname(1:startIndex(1)) 'MIPidx_' fname(startIndex(1)+1:end)];

            % if the channels were split in export, put in the channel index
            if ~isempty(regexp(fname,'_w%.4d','once'))
                fname = sprintf(fname, channel-1);
                
            % otherwise they were still split by makeMIP so we should
            % modify the filename
            else
                fname = [fname(1:end-4) sprintf('_w%.4d',channel-1) fname(end-3:end)];
            end
            
            fname = fullfile(dataDir, fname);
 
            % when I make MIPs I replace m or f by p, I want all positions
            % labeled the same and thats how Andor should have done it
            montageIdx = regexp(fname,'(m|f)[0-9]+','once');
            if ~isempty(montageIdx)
                fname(montageIdx) = 'p';
            end
            
            if ~exist(fname,'file')
                error(['MIPidx file does not exist: ' fname]);
                MIPidx = [];
            else
                MIPidx = imread(fname,time);
            end
            
            %disp(['loaded MIPidx ' fname]);
        end
        
        % process data
        %---------------------------------

        function debugInfo = extractData(this, dataDir, nucChannel, opts)
            % populate the data array
            %
            % extractData(dataDir, nuclearChannel)
            % extractData(dataDir, nuclearChannel, options)
            %
            % dataDir:              data directory
            % nuclearChannel        channel of nuclear segmentation
            % options:              struct with fields
            %
            % -dataChannels         channels for which to extract data
            %                       default is all
            % -tMax                 maximal time, default nTime
            % -------------------------------------------------------------
            % -cleanupOptions       options for cleanup
            % -nucShrinkage         number of pixels to shrink nuclear mask
            %                       by to get more accurate nuclear
            %                       readouts
            % -cytoplasmicLevels    extract approximate cytoplasmic levels
            % -cytoSize             width of cytoplasmic annulus
            % -fgChannel            channel whose inverse is used as a
            %                       background mask
            % -bgMargin             erosion size of the inverse fg 
            % -imopenBGsub          open size for imopenBGsub (if using it)
            % -------------------------------------------------------------
            % -MIPidxDir            directory of MIPidx files if desired
            % -segmentationDir      directory of segmentation, default
            %                       dataDir
            % -------------------------------------------------------------
            % -nuclearSegmentation  binary stack with 3rd dim time for
            %                       providing a segmentation by hand

            if nargin < 4
                opts = struct();
            end
            if ~isfield(opts,'dataChannels')
                opts.dataChannels = 1:this.nChannels;
            end
            if ~isfield(opts, 'tMax')
                opts.tMax = this.nTime;
            end
            % ----------------------------------------
            if ~isfield(opts, 'cleanupOptions')
                opts.cleanupOptions = struct('separateFused', true,...
                    'clearBorder',true, 'minAreaStd', 1, 'minSolidity',0);
            end
            if ~isfield(opts, 'nucShrinkage')
                opts.nucShrinkage = 0;
            end
            if ~isfield(opts,'cytoplasmicLevels')
                opts.cytoplasmicLevels = false;
            end
            if ~isfield(opts, 'cytoSize')
                opts.cytoSize = 10;
            end
            if ~isfield(opts, 'bgMargin')
                opts.bgMargin = 2;
            end
            if ~isfield(opts, 'imopenBGsub')
                opts.imopenBGsub = [];
            end
            % ----------------------------------------
            if ~isfield(opts, 'segmentationDir')
                opts.segmentationDir = dataDir;
            end
            if ~isfield(opts, 'MIPidxDir')
                opts.MIPidxDir = [];
            end
            
            this.dataChannels = opts.dataChannels;
            allChannels = unique([nucChannel opts.dataChannels]);

            % for the purpose of the current background subtraction
            if isfield(opts, 'fgChannel');
                allChannels = unique([allChannels opts.fgChannel]);
            else
                bgmask = [];
                fgmask = [];
            end
            
            % load all available segmentations
            seg = {};
            for ci = allChannels
                seg{ci} = this.loadSegmentation(opts.segmentationDir, ci);
            end
            
            % pass nuclear segmentation manually or load Ilastik standard
            if isfield(opts, 'nuclearSegmentation')            
                seg{nucChannel} = opts.nuclearSegmentation;
                if size(seg{nucChannel},3) ~= this.nTime
                    error('slice number of nuclear segmentation doesnt match number of time points');
                end
            end
            if isempty(seg{nucChannel})
                error('nuclear segmentation missing');
            end

            for ti = 1:opts.tMax

                % progress indicator
                fprintf('.');
                if mod(ti,60)==0
                    fprintf('\n');
                end
                
                % get background mask
                %--------------------------
                if isfield(opts, 'fgChannel');
                    fgmask = imclose(seg{opts.fgChannel}(:,:,ti),strel('disk',10));
                    bgmask = imerode(~fgmask,strel('disk', opts.bgMargin));
                    fgmask = imerode(fgmask,strel('disk', opts.bgMargin));
                else
                    warning('not using background mask');
                end
                
                % make clean nuclear mask
                %--------------------------
            
                nucmaskraw = seg{nucChannel}(:,:,ti);
                %nucmaskraw(bgmask) = false;
                if ~isfield(opts, 'nuclearSegmentation')  
                    nucmask = nuclearCleanup(nucmaskraw, opts.cleanupOptions);
                else
                    nucmask = nucmaskraw;
                end
                
                % make cytoplasmic mask 
                %----------------------
% OLD                
%                 if opts.cytoplasmicLevels
% 
%                     % watershedding inside the dilation
%                     dilated = imdilate(nucmask, strel('disk',opts.cytoSize));
%                     basin = imcomplement(bwdist(dilated));
%                     basin = imimposemin(basin, nucmask);
%                     L = watershed(basin);

%                     stats = regionprops(L, 'PixelIdxList');
%                     cytCC = struct('PixelIdxList', {cat(1,{stats.PixelIdxList})});
% 
%                     this.cellData(ti).cytLevel = zeros([numel(cytCC.PixelIdxList) numel(opts.dataChannels)]);
%                     this.cellData(ti).cytLevelAvg = zeros([1 numel(opts.dataChannels)]);
%                 end
                
                if opts.cytoplasmicLevels

                    % watershedding inside the dilation
                    dilated = imdilate(nucmask, strel('disk',opts.cytoSize));
                    basin = imcomplement(bwdist(dilated));
                    basin = imimposemin(basin, nucmask);
                    L = watershed(basin);
                    
                    % OLD: can change number of CC and mess things up:
                    % L(~dilated) = 0;    % exclude outside dilated nuclei
                    % L(nucmaskraw) = 0;  % exclude nuclei before cleanup
                    % L(~fgmask) = 0;      % exclude background 
                    
                    stats = regionprops(L, 'PixelIdxList');
                    cytCC = struct('PixelIdxList', {cat(1,{stats.PixelIdxList})});
                    
                    % replacement of exclusions above
                    for cci = 1:numel(cytCC.PixelIdxList)
                        CCPIL = cytCC.PixelIdxList{cci};
                        CCPIL = CCPIL(dilated(CCPIL));    % exclude outside dilated nuclei
                        CCPIL = CCPIL(~nucmaskraw(CCPIL));% exclude nuclei before cleanup
                        if ~isempty(fgmask)
                            cytCC.PixelIdxList{cci} = CCPIL(fgmask(CCPIL));% exclude background 
                        end
                    end
                    
                    this.cellData(ti).cytLevel = zeros([numel(cytCC.PixelIdxList) numel(opts.dataChannels)]);
                    this.cellData(ti).cytLevelAvg = zeros([1 numel(opts.dataChannels)]);
                else
                    cytCC = {};
                end
                
                % regionprops and indices from nuclear mask
                % follows cyt so we can shrink without introducing new var
                %-------------------------------------------
                
                % shrink instead of erode to preserve number of connected
                % components to not introduce mismatch between nuclear and
                % cytoplasmic mask
                if opts.nucShrinkage > 0
                    nucmask = bwmorph(nucmask,'shrink',opts.nucShrinkage);
                end
                
                % initialize cellData
                %---------------------------------------------------
                
                nucCC = bwconncomp(nucmask);
                stats = regionprops(nucCC, 'Area', 'Centroid');
                centroids = cat(1,stats.Centroid);
                areas = cat(1,stats.Area);
                nCells = nucCC.NumObjects;

                this.ncells(ti) = nCells;
                this.cellData(ti).XY = centroids;
                this.cellData(ti).area = areas;
                this.cellData(ti).nucLevel = zeros([nCells numel(opts.dataChannels)]);
                this.cellData(ti).nucLevelAvg = zeros([1 numel(opts.dataChannels)]);
                this.cellData(ti).background = zeros([1 numel(opts.dataChannels)]);
                
                % if is no MIPidx, analyze the MIP
                % taking the MIP of a single z-slice costs no time so
                % this works for single images
                if ~isempty(opts.MIPidxDir)
                    MIPidx = this.loadMIPidx(opts.MIPidxDir, nucChannel, ti);
                else
                    MIPidx = [];
                end
                % ti == 1 so it doesnt say it a hundred times
                if ti == 1 && isempty(MIPidx) 
                   warning('------------ NO MIPidx FOUND ------------');
                end
                
                % for background subtraction, median z-plane
                % this ifempty(MIPidx) is just so that it proceeds without
                % MIPidx if the MIPidx cannot be loaded
                if ~isempty(MIPidx)
                    zmed = median(MIPidx(nucmask));
                else
                    % if there is no MIPidx, it will use the MIP
                    % read out in its only plane
                    zmed = 1; 
                end
                
                % read out nuclear and cytoplasmic levels 
                %-----------------------------------------
                % zmed == 0 means no nuclei in mask so no cells in MIPidx
                % case
                if nCells > 0  % zmed > 0 
                    
                for cii = 1:numel(opts.dataChannels)
                    
                    %disp(['loading channel ' num2str(opts.dataChannels(cii))]);
                    imc = this.loadImage(dataDir, opts.dataChannels(cii), ti);
                    %disp(['size: ' num2str(size(imc))]);
                    
                    if size(nucmask) ~= [size(imc,1) size(imc,2)]
                        error(['nucmask size ' num2str(size(nucmask)) ' does not match image size ' num2str(size(imc))]);
                    end
                    
                    % if no MIPidx, just use MIP
                    if isempty(MIPidx)
                        imc = max(imc,[],3);
                    end

                    if ~isempty(opts.imopenBGsub)
                        disp('imopen bg sub');
                        imc = imc - imopen(imc,strel('disk',opts.imopenBGsub));
                    end
                    
                    % current low-tech background subtraction:
                    %-------------------------------------------------
                    % mean value of segmented empty space in the image
                    % or otherwise just min of image
                    if ~isempty(bgmask)
                        % disp('bg seg bg sub');
                        % if the background area is too small to be
                        % reliable (I'm saying < ~1% of total area), then
                        % just go with the previous value
                        % zmed = 0 means no more nuclei left
                        if sum(bgmask(:)) > 10^4
                            imcZmed = imc(:,:,zmed);
                            this.cellData(ti).background(cii) = mean(imcZmed(bgmask));
                        else
                            if ti > 1
                                this.cellData(ti).background(cii) = this.cellData(ti-1).background(cii);
                            else
                                this.cellData(ti).background(cii) = min(min(imc(:,:,zmed)));
                            end
                        end
                    else
                        this.cellData(ti).background(cii) = min(min(imc(:,:,zmed)));
                    end
                    
                    for cellidx = 1:nCells
                        
                        if ~isempty(MIPidx)
                            zi = median(MIPidx(nucCC.PixelIdxList{cellidx}));
                        else
                            zi = 1;
                        end
                    
                        nucPixIdx = nucCC.PixelIdxList{cellidx};
                        nucPixIdx = double(zi-1)*size(imc,1)*size(imc,2) + nucPixIdx;
                        this.cellData(ti).nucLevel(cellidx, cii) = mean(imc(nucPixIdx));
                        
                        if opts.cytoplasmicLevels
                            cytPixIdx = cytCC.PixelIdxList{cellidx};
                            cytPixIdx = double(zi-1)*size(imc,1)*size(imc,2) + cytPixIdx;
                            this.cellData(ti).cytLevel(cellidx, cii) = mean(imc(cytPixIdx));
                        end
                    end
                    
                    nL = this.cellData(ti).nucLevel(:,cii);
                    A = this.cellData(ti).area;
                    if opts.cytoplasmicLevels
                        cL = this.cellData(ti).cytLevel(:,cii);
                        
                        % calculate individual nuc/cyt ratio to exclude garbage
                        % (e.g. really bright spot on top of nuc can give
                        % ratios >> 2, which cannot be real
                        bg = this.cellData(ti).background(cii);
                        NCR = (nL-bg)./(cL-bg); % nuc/cyt ratio
                        
                        % cytoplasmic mask can be empty if a cell was tiny
                        % or junk was defined as a nucleus
                        idx = NCR < 3 & ~isnan(cL);
                        this.cellData(ti).cytLevelAvg(cii) = mean(cL(idx).*A(idx))/mean(A(idx));
                    else
                        idx = 1:numel(A);
                    end
                    this.cellData(ti).nucLevelAvg(cii) = mean(nL(idx).*A(idx))/mean(A(idx));
                end
                else
                    warning(['------------ NO CELLS AT T = ' num2str(ti) '------------']);
                end
            end
            fprintf('\n');
            
            debugInfo = struct('bgmask', bgmask, 'nucmask', nucmask, 'cytCC', cytCC);
        end
        
        function makeTimeTraces(this)
            % make time traces of levels in cellData
            %
            % makeTimeTraces()
            % 
            % populates Position.timeTraces

            nucLevelAvg = zeros([this.nTime numel(this.dataChannels)]);
            cytLevelAvg = zeros([this.nTime numel(this.dataChannels)]);

            nucLevelMed = zeros([this.nTime numel(this.dataChannels)]);
            cytLevelMed = zeros([this.nTime numel(this.dataChannels)]);

            bg = zeros([this.nTime numel(this.dataChannels)]);

            for ti = 1:numel(this.cellData)

                nucLevel = this.cellData(ti).nucLevel;
                cytLevel = this.cellData(ti).cytLevel;

                % averages computed in extractData excluding outliers
                nucLevelAvg(ti, :) = this.cellData(ti).nucLevelAvg;
                cytLevelAvg(ti, :) = this.cellData(ti).cytLevelAvg;
                
                idx = ~isnan(nucLevel) & ~isnan(cytLevel);
                A = this.cellData(ti).area;
                nucLevelMed(ti, :) = median(nucLevel(idx).*A(idx))/mean(A(idx));
                cytLevelMed(ti, :) = median(cytLevel(idx).*A(idx))/mean(A(idx));
                
                bg(ti) = this.cellData(ti).background;
            end
            
            this.timeTraces.nucLevelAvg = nucLevelAvg;
            this.timeTraces.cytLevelAvg = cytLevelAvg;
            
            this.timeTraces.nucLevelMed = nucLevelMed;
            this.timeTraces.cytLevelMed = cytLevelMed;
            
            this.timeTraces.background = bg;
        end
        
        % getter for dependent properties
        %---------------------------------
        
        function barefname = get.bareFilename(this)
            % filename without extension 

            barefname = this.filename(1:end-4);
        end
        
%         data was for compatibility with CellTracker, but is not used by
%         me or anyone else afaik and is causing problems
%
%         function data = get.data(this)
%             % output CellTracker style data array for first time point
%             % 
%             % data cols: x, y, area, -1, 
%             % then (nuclear, cytoplasmic) mean for each channel
%             
%             % this is just for old static stuff for now
%             % perhaps it can be removed alltogether
%             if this.nTime > 1
%                 data = [];
%                 return;
%             end
%             
%             cData = this.cellData(1);
%             if ~isempty(this.ncells)
%                 ncells = this.ncells(1);
%             else
%                 ncells = 0;
%             end
%             
%             if ~isfield(cData,'nucLevel') || isempty(cData.nucLevel)
%                 data = [];
%                 return;
%             end
%             
%             if isfield(cData,'cytLevel') && ~isempty(cData.cytLevel)
%                 cytLevel = cData.cytLevel;
%             else
%                 cytLevel = nan(size(cData.nucLevel));
%             end
% 
%             data = [cData.XY, cData.area,...
%                 -ones([ncells 1]), zeros([ncells 2*numel(this.dataChannels)])];
%             
%             for cii = 1:numel(this.dataChannels)
%                 data(:,3 + 2*cii) = cData.nucLevel(:,cii);%this.dataChannels(cii)
%                 data(:,4 + 2*cii) = cytLevel(:,cii);%this.dataChannels(cii)
%             end
%         end

        % setters
        %---------------------------------
        
        function setID(this, ID)
            
            this.ID = ID;
        end
    end
end