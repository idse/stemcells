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
        function this = Position(nChannels, filename)

            % matlab sucks
            if nargin == 0
                return
            end

            this.nChannels = nChannels;
            this.filename = filename;
            
            this.cellData = struct();
        end

        % saving and loading
        %---------------------------------

        function img = loadImage(this, dataDir, channels)
            % load position image
            %
            % img = loadImage(dataDir, channels)
            %
            % dataDir:  main data directory 
            % channels: desired channels to be loaded, leave empty for all
            %
            % img:      loaded image

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
            
            % Ilastik output data has xy transposed
            seg = permute(seg, [2 1 3]);
        end
        
        % process data
        %---------------------------------

        function extractData(this, dataDir, nuclearChannel, dataChannels)
            % populate the data array
            %
            % extractData(dataDir, nuclearChannel)
            % extractData(dataDir, nuclearChannel, dataChannels)
            %
            % dataDir:          data directory
            % nuclearChannel:   channel of nuclear segmentation
            % dataChannels:     channels for which to extract data
            %                   default is all

            seg = this.loadSegmentation(dataDir, nuclearChannel);
            
            if nargin < 4
                dataChannels = 1:this.nChannels;
            end
            this.dataChannels = dataChannels;
            
            % most time consuming: nuclearCleanup
            options = struct('separateFused', true, 'clearBorder',true,...
                    'minAreaStd', 1, 'minSolidity',0);
                
            nucmaskraw = seg;
            nucmask = nuclearCleanup(nucmaskraw, options);
            %cytmask = imdilate(nucmask, strel('disk',10));
            %cytmask(nucmaskraw)=0;
                
            CC = bwconncomp(seg);
            stats = regionprops(CC, 'Area', 'Centroid');
            centroids = cat(1,stats.Centroid);
            areas = cat(1,stats.Area);
            nCells = CC.NumObjects;

            this.ncells = nCells;
            this.cellData.XY = centroids;
            this.cellData.area = areas;
            this.cellData.nucLevel = zeros([nCells this.nChannels]);
            
            for cii = 1:numel(dataChannels)
                imc = this.loadImage(dataDir, dataChannels(cii));
                imc = imc - min(imc(:));
                for cellidx = 1:nCells
                    this.cellData.nucLevel(cellidx, cii) = mean(imc(CC.PixelIdxList{cellidx}));
                end
            end
        end
        
        % getter for dependent properties
        %---------------------------------
        
        function barefname = get.bareFilename(this)
            % filename without extension 

            barefname = this.filename(1:end-4);
        end
        
        function data = get.data(this)
            % output CellTracker style data array
            % 
            % data cols: x, y, area, -1, 
            % then (nuclear, cytoplasmic) mean for each channel
            
            if ~isfield(this.cellData,'nucLevel') || isempty(this.cellData.nucLevel)
                data = [];
                return;
            end
            
            if isfield(this.cellData,'cytLevel') && ~isempty(this.cellData.cytLevel)
                cytLevel = this.cellData.cytLevel;
            else
                cytLevel = nan(size(this.cellData.nucLevel));
            end

            data = [this.cellData.XY, this.cellData.area,...
                -ones([this.ncells 1]), zeros([this.ncells 2*this.nChannels])];
            
            for ci = 1:this.nChannels
                data(:,3 + 2*ci) = this.cellData.nucLevel(:,ci);
                data(:,4 + 2*ci) = cytLevel(:,ci);
            end
        end

        % setters
        %---------------------------------
        
        function setID(this, ID)
            
            this.ID = ID;
        end
    end
end