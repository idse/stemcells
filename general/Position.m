classdef Position < handle
    % Data class to store cell data in a field of view

    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------

    properties

        ncells              % number of cells
        density             % cell density

        cellData            % structure:
                            % - XY
                            % - area
                            % - nucLevel:   nCell x nChannels array
                            % - cytLevel
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
        data                % cell by cell data for colony
                            % data cols: x, y, area, -1, 
                            % (nuclear, cytoplasmic) mean for each channel
    end

    methods

        % constructor
        function this = Position(nChannels, filename, nTime)

            % matlab sucks
            if nargin == 0
                return
            end
            if ~exist('nTime','var')
                nTime = 1;
            end

            this.nChannels = nChannels;
            this.nTime = nTime;
            this.filename = filename;
            
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

            if exist('time','var') && time > 1
                error('todo : include reading for dynamic not Andor datasets');
            end
            
            if ~exist('channels','var')
                channels = 1:this.nChannels;
            end
            fname = fullfile(dataDir, this.filename);
            [~,~,ext] = fileparts(this.filename);
            
            if strcmp(ext,'.tif')

                info = imfinfo(fname);
                w = info.Width;
                h = info.Height;

                img = zeros([h w numel(channels)],'uint16');
                for cii = 1:numel(channels)
                    img(:,:,cii) = imread(fname, channels(cii));
                end
            else
                img = readStack(fname);
                img = img(:,:,channels); % yes this is very inefficient
                                         % I will fix it up later 
            end
        end

        function seg = loadSegmentation(this, dataDir, channel)
            % load segmentation
            %
            % seg = loadSegmentation(dataDir)
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

            s = strsplit(this.filename,'_.%\.4d|\.','DelimiterType','RegularExpression');            
            barefname = s{1};
            
             % for dynamic data the raw data will usually be a MIP
            % this code tries both
            listing = dir(fullfile(dataDir,[barefname '_*h5']));
            if isempty(listing)
                s = strsplit(this.filename,'_');
                barefname = sprintf([s{1} '_MIP_p%.4d'], this.ID);
                listing = dir(fullfile(dataDir,[barefname '_*h5']));
            end
            if isempty(listing)
                error(['segmentation not found in ' dataDir]);
            end
            
            seg = [];
            
            for i = 1:numel(listing)
                if ~isempty(strfind(listing(i).name,sprintf('_w%.4d',channel-1)))... 
                        || ~isempty(strfind(listing(i).name,sprintf('_c%d',channel)))
                    fname = fullfile(dataDir,listing(i).name);
                    seg = (squeeze(h5read(fname, '/exported_data')) == 2);
                end
            end
            
            if isempty(seg)
                error(['segmentation not found in ' dataDir, ', may be naming convention problem']);
            end
            
            % Ilastik output data has xy transposed
            seg = permute(seg, [2 1 3]);
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
            
            [startIndex,endIndex]  = regexp(this.filename,'_t([0-9]+|%.4d)');
            
            % remove time part of filename:
            % MIPidx is supposed to be single file for all times
            fname = this.filename([1:startIndex-1 endIndex+1:end]);
            startIndex = regexp(fname,'_');
            fname = [fname(1:startIndex(1)) 'MIPidx_' fname(startIndex(1)+1:end)];
            fname = fullfile(dataDir, fname);
            
            % if the channels are split in Andor format, put it in
            if ~isempty(regexp(fname,'_w%.4d','once'))
                fname = sprintf(fname, channel-1);
            end
            
            % when I make MIPs I replace m by p, I want all positions
            % labeled p and thats how Andor should have done it
            montageIdx = regexp(fname,'m[0-9]+','once');
            if ~isempty(montageIdx)
                fname(montageIdx) = 'p';
            end
            
            if ~exist(fname,'file')
                warning(['MIPidx file does not exist: ' fname]);
                MIPidx = [];
            else
                MIPidx = imread(fname,time);
            end
        end
        
        % process data
        %---------------------------------

        function extractData(this, dataDir, nuclearChannel, opts)
            % populate the data array
            %
            % extractData(dataDir, nuclearChannel)
            % extractData(dataDir, nuclearChannel, dataChannels)
            %
            % dataDir:              data directory
            % nuclearChannel:       channel of nuclear segmentation
            % options:              struct with fields
            % -dataChannels         channels for which to extract data
            %                       default is all
            % -cleanupOptions       options for cleanup
            % -cytoplasmicLevels    extract approximate cytoplasmic levels
            % -MIPidxDir            directory of MIPidx files if desired
            % -tMax                 maximal time, default nTime
            % -segmentationDir      directory of segmentation, default
            %                       dataDir

            if nargin < 4
                opts = struct();
            end
            if ~isfield(opts,'dataChannels')
                opts.dataChannels = 1:this.nChannels;
            end
            if ~isfield(opts, 'cleanupOptions')
                opts.cleanupOptions = struct('separateFused', true,...
                    'clearBorder',true, 'minAreaStd', 1, 'minSolidity',0);
            end
            if ~isfield(opts,'cytoplasmicLevels')
                opts.cytoplasmicLevels = false;
            end
            if ~isfield(opts, 'tMax')
                opts.tMax = this.nTime;
            end
            if ~isfield(opts, 'segmentationDir')
                opts.segmentationDir = dataDir;
            end
            if ~isfield(opts, 'MIPidxDir')
                MIPidxDir = [];
            end
            
            this.dataChannels = opts.dataChannels;
            
            seg = this.loadSegmentation(opts.segmentationDir, nuclearChannel);
            
            % make clean nuclear mask and initialize cellData
            %-----------------------------------------------------------

            for ti = 1:opts.tMax

                % progress indicator
                fprintf('.');
                if mod(ti,60)==0
                    fprintf('\n');
                end
                
                nucmaskraw = seg(:,:,1);
                nucmask = nuclearCleanup(nucmaskraw, opts.cleanupOptions);
                
                nucCC = bwconncomp(nucmask);
                stats = regionprops(nucCC, 'Area', 'Centroid');
                centroids = cat(1,stats.Centroid);
                areas = cat(1,stats.Area);
                nCells = nucCC.NumObjects;

                this.ncells(ti) = nCells;
                this.cellData(ti).XY = centroids;
                this.cellData(ti).area = areas;
                this.cellData(ti).nucLevel = zeros([nCells numel(opts.dataChannels)]);

                % make cytoplasmic mask 
                %----------------------

                if opts.cytoplasmicLevels

                    % watershedding inside the dilation
                    dilated = imdilate(nucmask, strel('disk',10));
                    basin = imcomplement(bwdist(dilated));
                    basin = imimposemin(basin, nucmask);
                    L = watershed(basin);
                    L(~dilated) = 0;
                    L(nucmask) = 0;
                    stats = regionprops(L, 'PixelIdxList');
                    cytCC = struct('PixelIdxList', {cat(1,{stats.PixelIdxList})});

                    this.cellData(ti).cytLevel = zeros([nCells numel(opts.dataChannels)]);
                end
                
                % read out nuclear and cytoplasmic levels 
                %-----------------------------------------

                for cii = 1:numel(opts.dataChannels)
                    
                    imc = this.loadImage(dataDir, opts.dataChannels(cii), ti);

                    % if is no MIPidx, analyze the MIP
                    % taking the MIP of a single z-slice costs no time so
                    % this works for single images
                    if ~isempty(MIPidxDir)
                        MIPidx = this.loadMIPidx(MIPidxDir, nuclearChannel, time);
                    else
                        MIPidx = [];
                    end
                    % this construction is just so that it proceeds without
                    % MIPidx if the MIPidx cannot be loaded
                    if isempty(MIPidx)
                        imc = max(imc,[],3);
                        imc = imc - min(imc(:));
                    end

                    for cellidx = 1:nCells
                        
                        if ~isempty(MIPidx)
                            zi = median(MIPidx(CC.PixelIdxList{cellidx}));
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
                end
            end
            fprintf('\n');
        end
        
        % getter for dependent properties
        %---------------------------------
        
        function barefname = get.bareFilename(this)
            % filename without extension 

            barefname = this.filename(1:end-4);
        end
        
        function data = get.data(this)
            % output CellTracker style data array for first time point
            % 
            % data cols: x, y, area, -1, 
            % then (nuclear, cytoplasmic) mean for each channel
            
            % this is just for old static stuff for now
            % perhaps it can be removed alltogether
            if this.nTime > 1
                data = [];
                return;
            end
            
            cData = this.cellData(1);
            ncells = this.ncells(1);
            
            if ~isfield(cData,'nucLevel') || isempty(cData.nucLevel)
                data = [];
                return;
            end
            
            if isfield(cData,'cytLevel') && ~isempty(cData.cytLevel)
                cytLevel = cData.cytLevel;
            else
                cytLevel = nan(size(cData.nucLevel));
            end

            data = [cData.XY, cData.area,...
                -ones([ncells 1]), zeros([ncells 2*numel(this.dataChannels)])];
            
            for cii = 1:numel(this.dataChannels)
                data(:,3 + 2*cii) = cData.nucLevel(:,this.dataChannels(cii));
                data(:,4 + 2*cii) = cytLevel(:,this.dataChannels(cii));
            end
        end

        % setters
        %---------------------------------
        
        function setID(this, ID)
            
            this.ID = ID;
        end
    end
end