classdef DynamicPositionAndor < Position
    % Data class to store cell data in a field of view

    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------

    properties

        nTime
        tPerFile
    end
    
    methods

        function this = DynamicPositionAndor(meta, ID)
            % constructor 
            %
            % DynamicPositionAndor(meta, ID)
            %
            % meta:     MetadataAndor object
            % ID:       position index
            
            % matlab sucks
            if nargin == 0
                return
            end
            
            this.nChannels = meta.nChannels;
            this.cellData = struct();
            
            this.nTime = meta.nTime;
            this.tPerFile = meta.tPerFile;
            this.ID = ID;
            
            % this is a clunky way to only put the actual position in the
            % filename format but leave the rest of the %.4d pieces
            filenameFormat = meta.filename;
            barefname = sprintf(filenameFormat,this.ID);
            barefname = barefname(1:end-2);
            filenameFormat(1:numel(barefname)) = barefname;
            this.filename = filenameFormat;
        end

        % saving and loading
        %---------------------------------

        function img = loadImage(this, dataDir, channel, time)
            % load image of colony
            %
            % img = loadImage(dataDir, channel, time)
            %
            % dataDir:  main data directory 
            % channel:  desired channel to be loaded
            %
            % img:      loaded image
            %
            % for now assume Andor format input, can be expanded later

            % fti : time index of file, e.g. if tPerFile = 2
            % Andor times start at 0, our time index and subti at 1
            % subti : time index within file
            %
            % so time=1: fti = 0, subti = 1
            % time=30: fti = 16, subti = 1
            
            warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning');
            
            fti = ceil(time/this.tPerFile) - 1;
            subti = this.tPerFile-rem(time,this.tPerFile);
            fname = fullfile(dataDir, sprintf(this.filename,fti,channel-1));
            
            info = imfinfo(fname);
            w = info.Width;
            h = info.Height;
            nZslices = numel(info)/this.tPerFile;
            
            ioffset = 1 + (subti-1)*nZslices;
            img = zeros([h w nZslices],'uint16');
            for i = 1:nZslices
                img(:,:,i) = imread(fname, ioffset + i);
            end
        end
        
        function MIPidx = loadMIPidx(this, dataDir, channel, time)
            % load MIP index of position
            %
            % MIPidx = loadMIPidx(this, dataDir, channel)
            %
            % dataDir:  main data directory 
            % channels: desired channels to be loaded, leave empty for all
            %
            % img:      loaded image
            
            s = strsplit(this.filename,'_');
            fname = sprintf([s{1} '_MIPidx_p%.4d_w%.4d.tif'], this.ID, channel-1);
            fname = fullfile(dataDir,fname);

            if ~exist(fname,'file')
                warning(['MIPidx file does not exist: ' fname]);
                MIPidx = [];
            else
                MIPidx = imread(fname,time);
            end
        end
        
        % process data
        %---------------------------------

        function extractData(this, dataDir, MIPdir, nucChannel, dataChannels)
            % populate the data array
            %
            % extractData(dataDir, MIPdir, nuclearChannel, dataChannels)
            %
            % dataDir:  main data directory 
            % MIPdir:   has MIP, segmentation, MIPidx
            
            seg = this.loadSegmentation(MIPdir, nucChannel);

            if nargin == 5
                this.dataChannels = dataChannels;
            else
                this.dataChannels = 1:this.nChannels;
            end
            
            % most time consuming step: separateFused
            options = struct('separateFused', true, 'clearBorder',true,...
                    'minAreaStd', 1, 'minSolidity',0);

            for ti = 1

                nucmaskraw = seg(:,:,ti);
                nucmask = nuclearCleanup(nucmaskraw, options);
                cytmask = imdilate(nucmask, strel('disk',10));
                cytmask(nucmaskraw)=0;

                % CONTINUE HERE
                
                CC = bwconncomp(nucmask);
                stats = regionprops(CC, 'Area', 'Centroid');
                centroids = cat(1,stats.Centroid);
                areas = cat(1,stats.Area);
                nCells = CC.NumObjects;

                this.ncells(ti) = nCells;
                this.cellData(ti).XY = centroids;
                this.cellData(ti).area = areas;
                this.cellData(ti).nucLevel = zeros([nCells this.dataChannels]);

                MIPidx = this.loadMIPidx(MIPdir, nucChannel, ti);
                
                for cii = 1:numel(dataChannels)
                        
                    imc  = this.loadImage(dataDir, dataChannels(cii), ti);
                    
                    for cellidx = 1:nCells
                        
                        nucPixIdx = CC.PixelIdxList{cellidx};
                        
                        if isempty(MIPidx)
                            % analyze MIP
                            if size(imc,3) > 1
                                imc = max(imc,[],3);
                            end
                            imc = imc - min(imc(:));
                        else
                            zi = median(MIPidx(CC.PixelIdxList{cellidx}));
                            nucPixIdx = double(zi-1)*size(imc,1)*size(imc,2) + nucPixIdx;
                        end

                        this.cellData(ti).nucLevel(cellidx, cii) = mean(imc(nucPixIdx));
                    end
                end
            end
        end
        
        % setters
        %---------------------------------
        
        function setID(this, ID)
            
            this.ID = ID;
        end
    end
end