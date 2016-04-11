classdef position < handle
    % Data class to store cell data in a field of view
    
    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------

    properties
        
        data            % cell by cell data for colony
                        % data cols: x, y, area, -1, 
                        % then (nuclear, cytoplasmic) mean for each channel

        ncells          % number of cells
        density         % cell density
        
        % different from cellTracker: 
        %------------------------------
        
        nChannels       % number of colors imaged
        filename   
    end
    
    properties (SetAccess = protected)
        ID              % identifier of colony 
    end
    
    properties (Dependent)
        
        bareFilename;   % filename without extension
    end

    methods
        
        % constructor
        function this = position(nChannels, filename)

            % matlab sucks
            if nargin == 0
                return
            end
            
            this.nChannels = nChannels;
            this.filename = filename;
        end

        % saving and loading
        %---------------------------------
        
        function img = loadImage(this, dataDir, channels)
            % load image of colony
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

        function seg = loadSegmentation(this, dataDir)
            % load colony segmentation
            %
            % seg = loadSegmentation(dataDir)
            %
            % dataDir:  main data directory (segmentation in subdir 'colonies')
            %
            % seg:      binary image of nuclei
            %
            % segmentation is expected to have same name as colony DAPI tif
            % image but with .tif -> '_Simple Segmentation.h5'
            %
            % by convention use green for foreground in Ilastik
            
            listing = dir(fullfile(dataDir,[this.bareFilename '_*h5']));
            seg = {};
                
            for i = 1:numel(listing)
                fname = fullfile(dataDir,listing(i).name);
                seg{i} = (squeeze(h5read(fname, '/exported_data')) == 2)';
            end
            
            seg = cat(3,seg{:});
        end
        
        % process data
        %---------------------------------
        
        function extractData(this, dataDir)
            % populate the data array
            %
            % extractData(dataDir)
            %
            % dataDir:  main data directory (colonies in subdir 'colonies')
            % 
            % columns in this.data:
            % x, y, area, -1, then (nuclear, cytoplasmic) mean for each channel
            
            seg = this.loadSegmentation(dataDir);
            if size(seg,3) > 1
                seg = any(seg,3);
            end
            %seg = imerode(seg, strel('disk',2));

            % most time consuming:
            options = struct('minAreaStd', 1, 'minSolidity',0);
            seg = separateFusedNuclei(seg, options);

            CC = bwconncomp(seg);
            stats = regionprops(CC, 'Area', 'Centroid');
            centroids = cat(1,stats.Centroid);
            areas = cat(1,stats.Area);
            nCells = CC.NumObjects;

            this.ncells = nCells;
            this.data = [centroids, areas, -ones([nCells 1]),...
                                        zeros([nCells 2*this.nChannels])];
            colim = this.loadImage(dataDir);
            
            for ci = 1:this.nChannels
                colimc = colim(:,:,ci);
                colimc = colimc - min(colimc(:));
                for cellidx = 1:nCells
                    this.data(cellidx, 3 + 2*ci) = mean(colimc(CC.PixelIdxList{cellidx}));
                end
            end
        end
        
        % getter for dependent properties
        %---------------------------------
        
        function barefname = get.bareFilename(this)
            % filename without extension 

            [~,barefname] = fileparts(this.filename);
        end
        
        % setters
        %---------------------------------
        
        function setID(this, ID)
            
            this.ID = ID;
        end
    end
end