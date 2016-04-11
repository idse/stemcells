classdef colony < handle
    % Data class to store a stem cell colony
    
    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------
    
    properties
        
        data            % cell by cell data for colony
                        % data cols: x, y, area, -1, 
                        % then (nuclear, cytoplasmic) mean for each channel

        ncells          % number of cells
        center          % x-y pixel coordinates of center (relative to btf)
        density         % cell density
        
        % different from cellTracker: 
        %------------------------------
        
        ID              % identifier of colony
        nChannels       % number of colors imaged
        radiusPixel     % radius of colony
        radiusMicron  
        radialAvg       % radial average 
        radialAvgSeg    % radial average based on segmentation
        radialStd       % radial standard deviation
        radialStdSeg  
        radialBinEdges  % radial bin edges in pixels
                        % assumed to be the same for segmented and
                        % unsegmented
        boundingBox     % [xmin xmax ymin ymax] into btf
    end
    
    properties (Dependent)
        radius;         % =radiusPixel for compatibility with CellTracker
        bareFilename;   % filename with variable extensions
    end
    
    methods
        
        % constructor
        function this = colony(nChannels, center, radiusPixel, radiusMicron, boundingBox)

            % matlab sucks
            if nargin == 0
                return
            end
            
            this.nChannels = nChannels;
            this.center = center;
            this.radiusPixel = radiusPixel;
            this.radiusMicron = radiusMicron;
            this.boundingBox = boundingBox;
        end
        
        % saving and loading
        %---------------------------------
        
        function img = loadImage(this, dataDir, channels)
            % load image of colony
            %
            % img = loadImage(dataDir, channels)
            %
            % dataDir:  main data directory (colonies in subdir 'colonies')
            % channels: desired channels to be loaded, leave empty for all
            %
            % img:      loaded image
            
            if ~exist('channels','var')
                channels = 1:this.nChannels;
            end
            
            fname = fullfile(dataDir, 'colonies', [this.bareFilename '.tif']);
            w = this.boundingBox(2)-this.boundingBox(1) + 1;
            h = this.boundingBox(4)-this.boundingBox(3) + 1;
            
            img = zeros([h w numel(channels)],'uint16');
            for cii = 1:numel(channels)
                img(:,:,cii) = imread(fname, channels(cii));
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

            fname = fullfile(dataDir,'colonies',[this.bareFilename '_c1_Simple Segmentation.h5']);
            seg = (squeeze(h5read(fname, '/exported_data')) == 2)';
        end
        
        function saveImage(this, img, dataDir, channels)
            % save image of colony
            %
            % img = saveImage(img, dataDir, channels)
            %
            % img:      colony image
            % dataDir:  main data directory (colonies in subdir 'colonies')
            % channels: desired channels to be saved, leave empty for all
            %           if specifying multiple channels, separate files
            %           will be made for each, if leaving empty a single
            %           stack will be saved with all channels

            if exist('channels', 'var')
                fname = fullfile(dataDir, 'colonies', [this.bareFilename '_c%d.tif']);
                for ci = channels
                    imwrite(img(:,:,ci), sprintf(fname, ci));
                end
            else
                fname = fullfile(dataDir, 'colonies', [this.bareFilename '.tif']);
                imwrite(img(:,:,1), fname);
                for ci = 2:this.nChannels
                    imwrite(img(:,:,ci), fname,'WriteMode','Append');
                end
            end
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
            this.data = [centroids, areas, -ones([nCells 1]), zeros([nCells 2*this.nChannels])];
            for ci = 1:this.nChannels
                colim = this.loadImage(dataDir, ci);
                colim = colim - min(colim(:));
                for cellidx = 1:nCells
                    this.data(cellidx, 3 + 2*ci) = mean(colim(CC.PixelIdxList{cellidx}));
                end
            end
        end
        
        function makeRadialAvgSeg(this)
            % create radial profile of segmented single cell data 
            % 
            % makeRadialAvgSeg(edges)
            %
            % edges:    edges for radial histogram in pixels
            %
            % populates this.radialAvgSeg where rows correspond to
            % radial bins, and columns to channels
            
            XY = this.data(:,1:2);
            XY(:,1) = XY(:,1) - mean(XY(:,1));
            XY(:,2) = XY(:,2) - mean(XY(:,2));
            r = sqrt(sum(XY(:,1:2).^2,2));
            [n,bini] = histc(r, this.radialBinEdges);

            nBins = numel(n)-1; % last bin we ignore (see doc histc)
            this.radialAvgSeg = zeros([nBins this.nChannels]);
            
            for ci = 1:this.nChannels    
                for i = 1:nBins
                    bindata = this.data(bini == i, 3 + 2*ci);
                    this.radialAvgSeg(i,ci) = mean(bindata);
                    this.radialStdSeg(i,ci) = std(bindata);
                end
            end
        end
        
        % getter for dependent properties
        %---------------------------------
        
        function radius = get.radius(this)
            % for compatibility with the CellTracker.colony
            radius = this.radiusPixel;
        end
        
        function barefname = get.bareFilename(this)
            % filename without extension 
            diam = this.radiusMicron*2;
            barefname = ['col_d' num2str(diam) '_id' num2str(this.ID)];
        end
    end
end