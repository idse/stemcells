classdef Colony < Position
    % Data class to store a stem cell colony
    
    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------
    
    properties

        center          % x-y pixel coordinates of center (relative to btf)
        well            % well located in
        
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
        
        %CM              % center of mass -> should use '.center'
        I               % second moments (inertia tensor)
    end
    
    properties (Dependent)
        radius;         % =radiusPixel for compatibility with CellTracker
    end
    
    methods
        
        % constructor
        function this = Colony(nChannels, center, radiusPixel, radiusMicron, boundingBox, well, nTime)
            % Colony(nChannels, center, radiusPixel, radiusMicron, boundingBox, well, nTime)

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

        function saveImage(this, img, dataDir, channels, MIP)
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

            if ~exist('MIP', 'var')
                MIP = true;
            end
            if MIP
                img = max(img,[],4);
            end

            if ~exist('channels', 'var') || isempty(channels)
                fname = fullfile(dataDir, [this.bareFilename '.tif']);
                imwrite(img(:,:,1), fname);
                for ci = 2:size(img, 3)
                    imwrite(img(:,:,ci),fname,'WriteMode','Append');
                    if ~MIP
                        error('all colors and no MIP is useless?');
                    end
                end
            else
                fname = fullfile(dataDir, [this.bareFilename '_c%d.tif']);
                for ci = channels
                    imwrite(img(:,:,ci,1), sprintf(fname, ci));
                    if ~MIP
                        for zi = 2:size(img, 4)
                            imwrite(img(:,:,ci,zi), sprintf(fname, ci),'WriteMode','Append');
                        end
                    end
                end
            end
        end

        % process data
        %---------------------------------

        function makeRadialAvgSeg(this, channels, normChannel, badidx)
            % create radial profile of segmented single cell data 
            % also calculate moment of inertia as measure of asymmetry
            % 
            % makeRadialAvgSeg()
            %
            % populates this.radialProfile.AvgSeg where rows correspond to
            % radial bins, and columns to channels
            
            if ~exist('channels','var') || isempty(channels)
                channels = 1:numel(this.dataChannels);
            end
            
            this.radialProfile = struct('BinEdges',[],'NucAvgSeg',[],...
                'NucStdSeg',[],'CytAvgSeg',[],'CytStdSeg',[],'NucCytRatio',[]);
            
            binWidthMicron = 10; % about two cell widths
            N = this.radiusMicron/binWidthMicron;
            binEdges = sqrt(linspace(0,this.radiusPixel^2,N+1));
            
            for ti = 1:numel(this.cellData)

                if ~exist('badidx','var') || isempty(badidx)
                    badidxti = false(size(this.cellData(ti).XY(:,1)));
                else
                    badidxti = badidx{ti};
                end
            
                if ~isempty(this.cellData(ti).XY)
    
                radialProf = struct('BinEdges',binEdges);
                
                XY = this.cellData(ti).XY;
                XY(:,1) = XY(:,1) - mean(XY(:,1));
                XY(:,2) = XY(:,2) - mean(XY(:,2));
                
                r = sqrt(sum(XY(:,1:2).^2,2));
                [n,bini] = histc(r, radialProf.BinEdges);
                nBins = numel(n)-1; % last bin we ignore (see doc histc)
                
                N = numel(this.dataChannels);
                radialProf.NucAvgSeg = zeros([nBins N]);
                radialProf.NucStdSeg = zeros([nBins N]);
                radialProf.CytAvgSeg = zeros([nBins N]);
                radialProf.CytStdSeg = zeros([nBins N]);
                radialProf.NucCytRatio = zeros([nBins N]);

                cyt = ~isempty(this.cellData(ti).cytLevel);
                
                for cii = channels
                    
                    for i = 1:nBins
                        
                        idx = bini == i & ~badidxti;
                        nucbindata = this.cellData(ti).nucLevel(idx, cii)...
                                            - this.cellData(ti).background(cii);

                        % make nuc:cyt before potentially normalizing nuc
                        % by DAPI levels below
                        if cyt 
                            cytbindata = this.cellData(ti).cytLevel(idx, cii)...
                                                - this.cellData(ti).background(cii);
                            radialProf.CytAvgSeg(i,cii) = mean(cytbindata);
                            radialProf.CytStdSeg(i,cii) = std(cytbindata);
                            
                            R = nucbindata./cytbindata;
                            %R(R < 0.3) = [];
                            %R(R > 2) = [];
                            radialProf.NucCytRatio(i,cii) = mean(R);
                        end
                        
                        % 180626 the behavior is different here from old code
                        % commented out below, where normChannel = 0
                        % normalized by cytoplasm, because here NucCytRatio
                        % is its own thing, this may affect Smad2 analysis
                        % if rerun
                        if exist('normChannel','var') && ~isempty(normChannel) && normChannel ~= 0
                            bg = this.cellData(ti).nucLevel(idx, normChannel)...
                                    - this.cellData(ti).background(normChannel);
                            nucbindata = nucbindata./bg;
                        end     
                        
                        radialProf.NucAvgSeg(i,cii) = mean(nucbindata);
                        radialProf.NucStdSeg(i,cii) = std(nucbindata);
                    end
                end
                this.radialProfile(ti) = radialProf;
            end
            end
            
            % make radial time traces
            for cii = channels
                
                this.timeTraces.radAvgNuc{cii} = zeros([this.nTime nBins]);
                this.timeTraces.radAvgCyt{cii} = zeros([this.nTime nBins]);
                this.timeTraces.radAvgRatio{cii} = zeros([this.nTime nBins]);

                this.timeTraces.moment{cii} = zeros([this.nTime 2]);
                
                for ti = 1:numel(this.cellData)
                    
                    % moment
                    XY = this.cellData(ti).XY - mean(this.cellData(ti).XY);
                    nucdata = this.cellData(ti).nucLevel(:,cii);
                    this.timeTraces.moment{cii}(ti,:) = sum(XY.*nucdata)/sum(nucdata);

                    % radial averages
                    this.timeTraces.radAvgNuc{cii}(ti,:) = this.radialProfile(ti).NucAvgSeg(:,cii);
                    this.timeTraces.radAvgCyt{cii}(ti,:) = this.radialProfile(ti).CytAvgSeg(:,cii);
                    this.timeTraces.radAvgRatio{cii}(ti,:) = this.radialProfile(ti).NucCytRatio(:,cii);
                end
            end
        end

%old version that was only for static data and no cytoplasmic levels
% USED FOR PAPER SMAD2 ANALYSIS
        function makeRadialAvgSegOld(this, channels, normChannel, badidx)
            % create radial profile of segmented single cell data 
            % 
            % makeRadialAvgSeg(normChannel, badidx)
            %
            % populates this.radialProfile.AvgSeg where rows correspond to
            % radial bins, and columns to channels
            
            nucLevel = this.cellData.nucLevel;
            if exist('normChannel','var')
                if normChannel ~= 0
                    nucLevel = nucLevel./this.cellData.nucLevel(:,normChannel);
                else
                    nucLevel = nucLevel./this.cellData.cytLevel;
                end
            end
            if exist('badidx','var')
                nucLevel = nucLevel(~badidx,:);
            end
            if ~exist('channels','var')
                channels = 1:numel(this.dataChannels);
            end
            
            XY = this.cellData.XY(~badidx,:);
            XY(:,1) = XY(:,1) - mean(XY(:,1));
            XY(:,2) = XY(:,2) - mean(XY(:,2));
            r = sqrt(sum(XY(:,1:2).^2,2));
            
            if isfield(this.radialProfile, 'BinEdges')
                binEdges = this.radialProfile.BinEdges;
            else
                binWidthMicron = 10; % about two cell widths
                N = this.radiusMicron/binWidthMicron;
                binEdges = sqrt(linspace(0,this.radiusPixel^2,N+1));
            end
            
            [n,bini] = histc(r, binEdges);
            nBins = numel(n)-1; % last bin we ignore (see doc histc)
            N = numel(channels);
            
            if ~isfield(this.radialProfile, 'NucAvgSeg')
                this.radialProfile.NucAvgSeg = zeros([nBins N]);
                this.radialProfile.NucStdSeg = zeros([nBins N]);
            end
            this.radialProfile.BinEdges = binEdges;
            
            for cii = channels
                for i = 1:nBins
                    bindata = nucLevel(bini == i, cii);
                    this.radialProfile.NucAvgSeg(i,cii) = mean(bindata);
                    this.radialProfile.NucStdSeg(i,cii) = std(bindata);
                end
            end
        end

        function makeRadialAvgNoSeg(this, colimg, colnucmask, colcytmask, colmargin, ti)
            
            % makeRadialAvgNoSeg(this, colimg, colnucmask, colcytmask, colmargin, ti)
            % colcytmask can be left empty
            
            if ~exist('ti','var')
                ti = 1;
            end
            
            % masks for radial averages
            [radialMaskStack, edges] = makeRadialBinningMasks(...
                this.radiusPixel, this.radiusMicron, colmargin);
            % colType is used when multiple size colonies are processed
            % together, not here
            colType = 1; 
            
            b = this.boundingBox(ti,:);
            % crop the maskstack if the colony is clipped 
            % right now only clipped left and top
            
            xmin = max(2-b(1),1);
            xmax = min(b(2), size(colimg,2)) + xmin - 1;

            ymin = max(2-b(3),1);
            ymax = min(b(4), size(colimg,1)) + ymin - 1;

            radialMaskStack{colType} = radialMaskStack{colType}(ymin:ymax,xmin:xmax,:);
%             crop = [max(2-b(1),1), min(b(2)-b(1)+1,size(radialMaskStack{colType},2)),...
%                     max(2-b(3),1), min(b(4)-b(3)+1,size(radialMaskStack{colType},1))];
%             radialMaskStack{colType} = radialMaskStack{colType}(crop(3):crop(4),crop(1):crop(2),:);
%             
            xres = this.radiusMicron / this.radiusPixel;
            if isempty(colcytmask)
                colcytmask = imdilate(colnucmask,strel('disk',round(5/xres)))-colnucmask;
            end
            
            % do radial binning
            N = size(radialMaskStack{colType},3);
            nucradavg = zeros([N this.nChannels]);
            nucradstd = zeros([N this.nChannels]);
            cytradavg = zeros([N this.nChannels]);
            cytradstd = zeros([N this.nChannels]);
            
            % store bin edges, to be reused by segmented profiles later
            this.radialProfile(ti).BinEdges = edges{colType};

            if this.nChannels > size(colimg,3) 
                nChannels = size(colimg,3);
                if ti == 1
                    warning('img is missing channel');
                end
            else
                nChannels = this.nChannels;
            end

            for ri = 1:N
                % for some reason linear indexing is faster than binary
                colnucbinmask = find(radialMaskStack{colType}(:,:,ri) & colnucmask);
                colcytbinmask = find(radialMaskStack{colType}(:,:,ri) & colcytmask);
                
                % we only want to measure if there is something there
                npix = sum(sum((colnucmask | colcytmask) & radialMaskStack{colType}(:,:,ri)));
                npixTot = sum(sum(radialMaskStack{colType}(:,:,ri)));
                
                if npix/npixTot > 0.1

                    for ci = 1:nChannels

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
%                 else
%                     fprintf('x');
                end
            end
            
            this.radialProfile(ti).NucAvg = nucradavg;
            this.radialProfile(ti).NucStd = nucradstd;
            this.radialProfile(ti).CytAvg = cytradavg;
            this.radialProfile(ti).CytStd = cytradstd;
        end

        function calculateMoments(this, colimg)%, colnucmask)
            
            siz = size(colimg);
            C = round(siz(1:2)/2);
            [X,Y] = meshgrid(1:siz(1),1:siz(2));
            X = X - C(1); Y = Y - C(2);

            this.CM = {};
            this.I = {};
            
            for ci = 1:this.nChannels
                
                % for testing
                %X0 = 200; Y0 = 300; R0 = 20;
                %M = (X-X0).^2 + (X-X0).*(Y-Y0) + (Y-Y0).^2 < R0^2;

                % background intensity outside colony
                R = round(this.radiusPixel);
                mask = X.^2 + Y.^2 < R.^2;
                
                M = colimg(:,:,ci);
                bg = mean(M(~mask));
                M = double(M).*mask;

                % weight for center of mass
                Mp = M - bg;
                W = double(Mp)./sum(Mp(:));
                this.CM{ci} = [sum(sum(W.*X)) sum(sum(W.*Y))];

                % inertia tensor
                Xp = X - this.CM{ci}(1);
                Yp = Y - this.CM{ci}(2);
                this.I{ci} = [   sum(sum(W.*Xp.*Xp)) sum(sum(W.*Yp.*Xp));
                        sum(sum(W.*Yp.*Xp)) sum(sum(W.*Yp.*Yp))];
                %[V,D] = eig(I);
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