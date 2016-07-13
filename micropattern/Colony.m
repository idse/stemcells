classdef Colony < Position
    % Data class to store a stem cell colony
    
    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------
    
    properties

        center          % x-y pixel coordinates of center (relative to btf)
        well % well located in
        % different from cellTracker: 
        %------------------------------
        
        radiusPixel     % radius of colony
        radiusMicron
        
        radialProfile   % struct with fields
                        % NucAvg    : radial average 
                        % NucAvgSeg : radial average based on segmentation
                        % NucStd    : radial standard deviation
                        % NucStdSeg  
                        % CytAvg    : cytoplasmic
                        % CytAvgSeg 
                        % CytStd
                        % CytStdSeg 
                        % BinEdges  : radial bin edges in pixels
                        % assumed to be the same for segmented and
                        % unsegmented
        boundingBox     % [xmin xmax ymin ymax] into btf
    end
    
    properties (Dependent)
        radius;         % =radiusPixel for compatibility with CellTracker
    end
    
    methods
        
        % constructor
        function this = Colony(nChannels, center, radiusPixel, radiusMicron, boundingBox, well, time)

            % matlab sucks
            if nargin == 0
                return
            end
            if ~exist('nTime','var')
                nTime = 1;
            end
            
            this.nTime = nTime;
            this.nChannels = nChannels;
            this.center = center;
            this.radiusPixel = radiusPixel;
            this.radiusMicron = radiusMicron;
            this.boundingBox = boundingBox;
            this.well = well; 
            d = num2str(2*radiusMicron);
            this.filename = ['col_d' num2str(d) '_id' num2str(this.ID) '.tif'];
        end
        
        % saving and loading
        %---------------------------------

        function saveImage(this, img, dataDir, channels)
            % save image of colony
            %
            % img = saveImage(img, dataDir, channels)
            %
            % img:      colony image
            % dataDir:  main data directory 
            % channels: desired channels to be saved, leave empty for all
            %           if specifying multiple channels, separate files
            %           will be made for each, if leaving empty a single
            %           stack will be saved with all channels

            if exist('channels', 'var')
                fname = fullfile(dataDir, [this.bareFilename '_c%d.tif']);
                for ci = channels
                    imwrite(img(:,:,ci), sprintf(fname, ci));
                end
            else
                fname = fullfile(dataDir, [this.bareFilename '.tif']);
                imwrite(img(:,:,1), fname);
                for ci = 2:this.nChannels
                    imwrite(img(:,:,ci), fname,'WriteMode','Append');
                end
            end
        end

        % process data
        %---------------------------------

        function makeRadialAvgSeg(this)
            % create radial profile of segmented single cell data 
            % 
            % makeRadialAvgSeg(edges)
            %
            % edges:    edges for radial histogram in pixels
            %
            % populates this.radialProfile.AvgSeg where rows correspond to
            % radial bins, and columns to channels
            
            XY = this.cellData.XY;
            XY(:,1) = XY(:,1) - mean(XY(:,1));
            XY(:,2) = XY(:,2) - mean(XY(:,2));
            r = sqrt(sum(XY(:,1:2).^2,2));
            [n,bini] = histc(r, this.radialProfile.BinEdges);

            nBins = numel(n)-1; % last bin we ignore (see doc histc)
            N = numel(this.dataChannels);
            this.radialProfile.NucAvgSeg = zeros([nBins N]);
            this.radialProfile.NucStdSeg = zeros([nBins N]);
            
            for cii = 1:N
                for i = 1:nBins
                    bindata = this.cellData.nucLevel(bini == i, cii);
                    this.radialProfile.NucAvgSeg(i,cii) = mean(bindata);
                    this.radialProfile.NucStdSeg(i,cii) = std(bindata);
                end
            end
        end

        % getter for dependent properties
        %---------------------------------
        
        function radius = get.radius(this)
            % for compatibility with the CellTracker.colony
            radius = this.radiusPixel;
        end
        
        % setters
        %---------------------------------
        
        function setID(this, ID)
            
            this.ID = ID;
            d = num2str(2*this.radiusMicron);
            this.filename = ['col_d' num2str(d) '_id' num2str(this.ID) '.tif'];
        end
    end
end