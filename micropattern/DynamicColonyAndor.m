classdef DynamicColonyAndor < DynamicPositionAndor
    % Data class to store cell data in a field of view

    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------

    properties

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
    end
    
    methods

        function this = DynamicColonyAndor(meta, ID, radiusMicron)
            % constructor 
            %
            % DynamicPositionAndor(meta, ID)
            %
            % meta:     MetadataAndor object
            % ID:       position index
            
            this = this@DynamicPositionAndor(meta, ID);
            
            % matlab sucks
            if nargin == 0
                return
            end
            
            this.radiusMicron = radiusMicron;
            this.radiusPixel = round(radiusMicron / meta.xres);
        end
        
        % process data
        %---------------------------------

        function makeRadialAvgSeg(this)
            % create radial profile of segmented single cell data 
            % 
            % makeRadialAvgSeg()
            %
            % populates this.radialProfile.AvgSeg where rows correspond to
            % radial bins, and columns to channels
            
            this.radialProfile = struct('BinEdges',[],'NucAvgSeg',[],...
                'NucStdSeg',[],'CytAvgSeg',[],'CytStdSeg',[],'NucCytRatio',[]);
            
            binWidthMicron = 10; % about two cell widths
            N = this.radiusMicron/binWidthMicron;
            binEdges = sqrt(linspace(0,this.radiusPixel^2,N+1));
            
            for ti = 1:numel(this.cellData)

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
                
                for cii = 1:N
                    for i = 1:nBins
                        
                        idx = bini == i;
                        nucbindata = this.cellData(ti).nucLevel(idx, cii)...
                                            - this.cellData(ti).background(cii);
                        radialProf.NucAvgSeg(i,cii) = mean(nucbindata);
                        radialProf.NucStdSeg(i,cii) = std(nucbindata);
                        
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
                    end
                end
                this.radialProfile(ti) = radialProf;
            end
            end
            
            % make radial time traces
            for cii = 1:N
                
                this.timeTraces.radAvgNuc{cii} = zeros([this.nTime nBins]);
                this.timeTraces.radAvgCyt{cii} = zeros([this.nTime nBins]);
                this.timeTraces.radAvgRatio{cii} = zeros([this.nTime nBins]);
                
                for ti = 1:numel(this.cellData)
                    this.timeTraces.radAvgNuc{cii}(ti,:) = this.radialProfile(ti).NucAvgSeg(:,cii);
                    this.timeTraces.radAvgCyt{cii}(ti,:) = this.radialProfile(ti).CytAvgSeg(:,cii);
                    this.timeTraces.radAvgRatio{cii}(ti,:) = this.radialProfile(ti).NucCytRatio(:,cii);
                end
            end
        end
    end
end