classdef cellStats < handle
	% Class to compare cell fate statistics

    % ---------------------
    % Idse Heemskerk, 2017
    % ---------------------

    properties
        
        % meta
        area
        conditions
        channelLabel
        markerChannels
        normalizeChannel
        
        % basic stats
        nucLevel 
        lim
        scale
        offset
        bins
        histograms

        % clusters
        clustermodel
        clusters
        nClusters
        clusterHistograms
        cutoff
    end
    
    methods

        % constructor
        function this = cellStats(allData, meta, markerChannels, normalizeChannel)
            
             if exist('normalizeChannel','var')
                this.normalizeChannel = normalizeChannel;
            else
                this.normalizeChannel = [];
            end

            this.conditions = meta.conditions;
            this.channelLabel = meta.channelLabel;
            this.markerChannels = markerChannels;
            
            % area in cm^2
            this.area = meta.xSize*meta.ySize*meta.xres*meta.yres/10^8;
            
            % make table of nuclear values
            for i = 1:numel(allData)
                
                if ~isempty(this.normalizeChannel)
                    N = allData{i}.cellData.nucLevel(:,normalizeChannel);
                    this.nucLevel{i} = bsxfun(@rdivide, allData{i}.cellData.nucLevel, N);
                else
                    this.nucLevel{i} = allData{i}.cellData.nucLevel; 
                end
            end
            
            % set intensity limits
            this.setLimits();
        end
        
        function setLimits(this)
            % determine intensity limits for distributions and plots
            
            nucLevelAll = [];
            for i = 1:numel(this.nucLevel)
                nucLevelAll = cat(1, nucLevelAll, this.nucLevel{i});
            end
    
            tol = 0.001;
            this.lim = {};
            this.scale = zeros([1 4]);
            this.offset = zeros([1 4]);
            for i = 1:size(nucLevelAll, 2)
                A = nucLevelAll(:,i);
                this.scale(i) = (max(A) - min(A));
                this.offset(i) = min(A);
                this.lim{i} = stretchlim(mat2gray(A),tol)*(max(A) - min(A)) + min(A);
                this.lim{i}(1) = min(A);
            end
            
%             for i = 1:numel(this.nucLevel)
%                 this.lim{i} = this.lim{i} - this.offset(i);
%                 for j = 1:size(nucLevelAll, 2)
%                     this.nucLevel{i}(:,j) = this.nucLevel{i}(:,j) - this.offset(j);
%                 end
%             end
        end

        function makeHistograms(this, nBins)
            % make histograms
            %
            % makeHistograms(nBins)
            
            this.histograms = {};
            this.bins = {};

            for fi = 1:numel(this.conditions)
                for channelIndex = 1:numel(this.lim)

                    this.bins{channelIndex} = linspace(this.lim{channelIndex}(1), this.lim{channelIndex}(2), nBins);
                    n = histc(this.nucLevel{fi}(:,channelIndex), this.bins{channelIndex});
                    this.histograms{fi, channelIndex} = n;
                end
            end
        end
        
        %-----------------------------------------------------------------
        % clusters
        %-----------------------------------------------------------------
        
        function [clusterLabels, nucInorm] = makeClusters(this, k, regularizationVal, cutoff)
            
            this.cutoff = cutoff;
            this.nClusters = k;
            
            mc = this.markerChannels;
            
            % combine nuclear values of all conditions
            nucI = this.nucLevel{1}(:,mc);
            condIdx = ones([size(nucI,1) 1],'uint8');
            for fi = 2:numel(this.conditions)
                nucI = cat(1, nucI, this.nucLevel{fi}(:,mc));
                condIdx = cat(1, condIdx, fi*ones([size(this.nucLevel{fi}(:,mc),1) 1],'uint8'));
            end

            % remove outliers and normalize the values
            nucInorm = nucI;
            for i = 1:numel(mc)
                outliers = nucI(:,i) > this.lim{mc(i)}(2);
                %nucI(outliers,:) = [];
                %nucInorm(outliers,:) = [];
                nucInorm(:,i) = nucInorm(:,i) - this.lim{mc(i)}(1);
                nucInorm(:,i) = nucInorm(:,i)./max(nucInorm(~outliers,i));
            end

            % fit Gaussian mixture model to all data
            this.clustermodel = fitgmdist(nucInorm(~outliers,:),k,'CovarianceType','full','SharedCovariance',false,...
                    'RegularizationValue', regularizationVal);
                
            % cluster all
            [clusterLabels, ~, clusterProb, ~] = cluster(this.clustermodel, nucInorm);
            clusterLabels(outliers) = 0;
            
            % cluster every individual conditions
            this.clusters = {};
            for fi = 1:numel(this.conditions)
                
                idx = condIdx == fi;
                this.clusters{fi} = struct('clusterLabels', clusterLabels(idx),...
                                            'clusterProb',  clusterProb(idx,:));

                % add mean and std
                for k = 1:this.nClusters
                    clusterIdx = this.getClusterIdx(fi, k);
                    this.clusters{fi}.mean{k} = mean(this.nucLevel{fi}(clusterIdx,:));
                    this.clusters{fi}.std{k} = std(this.nucLevel{fi}(clusterIdx,:));
                end
            end

%             % cluster every individual conditions, by applying model to
%             each case separately -> didn't work that well
%             this.clusters = {};
%             for fi = 1:numel(this.conditions)
% 
%                 Xc = this.nucLevel{fi}(:,this.markerChannels);
%                 Xcnorm = Xc;
%                 for i = 1:numel(this.markerChannels)
%                     Xcnorm(:,i) = Xcnorm(:,i) - this.lim{this.markerChannels(i)}(1);
%                     Xcnorm(:,i) = Xcnorm(:,i)./max(Xcnorm(:,i));
%                 end
% 
%                 [clusterLabels, ~, clusterProb, ~] = cluster(this.clustermodel, Xcnorm);
%                 this.clusters{fi} = struct('clusterLabels', clusterLabels, 'clusterProb',  clusterProb);
%                 
%                 % add mean and std
%                 for k = 1:this.nClusters
%                     clusterIdx = this.getClusterIdx(fi, k);
%                     this.clusters{fi}.mean{k} = mean(this.nucLevel{fi}(clusterIdx,:));
%                     this.clusters{fi}.std{k} = std(this.nucLevel{fi}(clusterIdx,:));
%                 end
%             end
        end
    
        function clusterIdx = getClusterIdx(this, conditionIdx, clusterLabel)
            
            C = this.clusters{conditionIdx};
            clusterIdx = C.clusterLabels == clusterLabel;
            goodIdx = max(C.clusterProb,[],2) > this.cutoff;
            clusterIdx = clusterIdx & goodIdx;
        end
        
        function makeClusterHistograms(this)
            
            for clusterLabel = 1:this.nClusters
                for fi = 1:numel(this.conditions)
                    
                    clusterIdx = this.getClusterIdx(fi, clusterLabel);
                    
                    for channelIndex = 1:numel(this.lim)

                        n = histc(this.nucLevel{fi}(clusterIdx,channelIndex), this.bins{channelIndex});
                        this.clusterHistograms{fi, channelIndex, clusterLabel} = n;
                    end
                end
            end
        end
        
        %-----------------------------------------------------------------
        % summary
        %-----------------------------------------------------------------

        function printSummary(this)

            clusterLabel = [];
            intensitySummary(clusterLabel);
            
            %-----------------------------------

            if ~isempty(this.clusters)
                
                for clusterLabel = 1:this.nClusters
                    intensitySummary(clusterLabel)
                end
            
                %-----------------------------------
                
                width = 12;
                conditionstr = this.conditions;
                for i = 1:numel(this.conditions)
                    N = numel(conditionstr{i}) + 6;
                    conditionstr{i} = pad(conditionstr{i}, N);
                end
                t = [pad('', width) conditionstr{:}];

                disp(' ');
                disp(t)
                disp(repmat('-',1,numel(t)));
                
                % cell numbers
                s = pad(['Ncells '], width);
                valstr = [];
                for fi = 1:numel(this.conditions)

                    N = numel(conditionstr{fi});
                    Nt = size(this.nucLevel{fi},1);
                    valstr = [valstr pad(sprintf('%d', Nt), N)];
                end
                disp([s valstr]);
                
                % density
                s = pad(['cells/cm^2'], width);
                valstr = [];
                for fi = 1:numel(this.conditions)

                    N = numel(conditionstr{fi});
                    rho = size(this.nucLevel{fi},1)/this.area;
                    valstr = [valstr pad(sprintf('%.1e', rho), N)];
                end
                disp([s valstr]);

                % fractions in cluster
                for clusterLabel = 1:this.nClusters
                    
                    s = pad(['cluster ' num2str(clusterLabel)], width);
                    valstr = [];
                    for fi = 1:numel(this.conditions)
                        
                        N = numel(conditionstr{fi});
                        Nt = size(this.clusters{fi}.clusterLabels,1);
                        Nc = sum(this.getClusterIdx(fi,clusterLabel));
                        fraction = Nc/Nt;
                        valstr = [valstr pad(sprintf('%.2f', fraction), N)];
                    end
                    disp([s valstr]);
                end
            end

            function intensitySummary(clusterLabel)

                width = 16;
                N = max(cellfun(@numel, this.conditions)) + 5;
                labels = pad(this.channelLabel(this.markerChannels), width);
                if isempty(clusterLabel)
                    t = [pad('',N) labels{:}];
                else
                    t = [pad(['cluster ' num2str(clusterLabel)], N) labels{:}];
                end
                
                disp(' ');
                disp(t)
                disp(repmat('-',1,numel(t)));
                for fj = 1:numel(this.conditions)

                    if isempty(clusterLabel)
                        idx = true([size(this.nucLevel{fj},1) 1]);
                    else
                        idx = this.getClusterIdx(fj, clusterLabel);
                    end
                    
                    s = pad(this.conditions{fj}, N);

                    meanval = mean(this.nucLevel{fj}(idx, this.markerChannels));
                    stdval = std(this.nucLevel{fj}(idx, this.markerChannels));

                    valstr = [];
                    for channel = 1:numel(this.markerChannels)
                        valstr = [valstr pad(sprintf('%.0f (%.0f)', meanval(channel), stdval(channel)), width)];
                    end
                    disp([s valstr])
                end
            end
        end
        
        %-----------------------------------------------------------------
        % visualization
        %-----------------------------------------------------------------
        
        function makeScatterPlot(this, conditionIdx, channelIdx, showClusters)
            
            if ~exist('showClusters','var')
                showClusters = false;
            end
            
            c1 = channelIdx(1);
            c2 = channelIdx(2);
            
            cdi = conditionIdx;
            
            hold on

            A = this.nucLevel{cdi}(:,c1);
            B = this.nucLevel{cdi}(:,c2);

            fs = 15;
            colormap(lines(this.clustermodel.NumComponents));

            if ~showClusters
                scatter(A, B, 1,'.');
                disp('no clusters');
            else
                clusterLabels = this.clusters{cdi}.clusterLabels;
                idx = max(this.clusters{cdi}.clusterProb,[],2) > this.cutoff;
                scatter(A, B, 1, '.','MarkerEdgeColor', 0.6*[1 1 1]);
                
                scatterlegend = {'unassigned'};
                
                for i = 1:max(clusterLabels)
                    idxx = clusterLabels == i & idx;
                    scatter(A(idxx), B(idxx), 1, clusterLabels(idxx),'.');
                    scatterlegend = [scatterlegend, ['cluster ' num2str(i)]];
                end
                
                legend(scatterlegend)
            end

            xlabel(this.channelLabel{c1}, 'FontSize',fs, 'FontWeight','Bold');
            ylabel(this.channelLabel{c2}, 'FontSize',fs, 'FontWeight','Bold');
            title(this.conditions{cdi},...
                'Interpreter','none', 'FontSize',fs, 'FontWeight','Bold');

            axis([this.lim{c1}' this.lim{c2}']);

            set(gcf,'color','w');
            set(gca, 'LineWidth', 2);
            set(gca,'FontSize', fs)
            set(gca,'FontWeight', 'bold')

            hold off
        end
        
        %-----------------------------------------------------------------
        
        function plotDistributionComparison(this, channelIndex, cumulative, clusterLabel)

            if ~exist('cumulative','var')
                cumulative = false;
            end
            if ~exist('clusterLabel','var')
                clusterLabel = [];
            end

            titlestr = this.channelLabel{channelIndex};
            
            if ~isempty(clusterLabel)
                nall = [this.clusterHistograms{:, channelIndex, clusterLabel}];
                titlestr = [titlestr ';  cluster ' num2str(clusterLabel)];
            else
                nall = [this.histograms{:, channelIndex}];
            end
            
            dist = bsxfun(@rdivide, nall, sum(nall,1));
            if cumulative
                dist = cumsum(dist);
            end

            [x,y] = histForBarlikePlot(this.bins{channelIndex}, dist);

            figure, 

            fs = 15;
            plot(x,y,'LineWidth',2);
            title(titlestr, 'FontSize',fs, 'FontWeight','Bold');

            xlabel('Intensity', 'FontSize',fs, 'FontWeight','Bold');
            ylabel('Frequency', 'FontSize',fs, 'FontWeight','Bold');

            set(gcf,'color','w');
            set(gca, 'LineWidth', 2);
            set(gca,'FontSize', fs)
            set(gca,'FontWeight', 'bold')

            xlim(this.lim{channelIndex});

            if cumulative
                legend(this.conditions, 'Location','SouthEast');
                ylim([0 1.1]);
            else
                legend(this.conditions, 'Location','NorthEast');
            end
        end
        
        %-----------------------------------------------------------------

        function showClusters(this, channelIdx)
            % plot Gaussian covariance as ellipse
            % UNFINISHED
            % NEED TO STORE VALUES TO SCALE BACK TO ORIGINAL COORDINATES

            hold on
            fs = 15;
            gmfit = this.clustermodel;

            c1 = channelIdx(1);
            c2 = channelIdx(2);
            
            ci1 = this.markerChannels == c1;
            ci2 = this.markerChannels == c2;
            
            %         % means of gaussians
            %         scatter(gmfit.mu(:,1), gmfit.mu(:,2),100,'.k')

            % covariance ellipse and cluster labels
            for i = 1:this.clustermodel.NumComponents

                text(gmfit.mu(i,ci1), gmfit.mu(i,ci2), num2str(i),'FontSize',fs, 'FontWeight','Bold');

                [V,D] = eig(gmfit.Sigma(:,:,i));
                t = linspace(0,2*pi);

                % major and minor axis length
                lmin = sqrt(D(1,1));
                lmaj = sqrt(D(2,2));

                majV = V(:,2);
                phi = atan(majV(2)/majV(1));
                R = [cos(phi) sin(phi); -sin(phi) cos(phi)];

                %quiver(gmfit.mu(i,1), gmfit.mu(i,2), majV(1), majV(2));

                X = lmaj*cos(t);
                Y = lmin*sin(t);

                Xp = R(1,1)*X + R(2,1)*Y + gmfit.mu(i,ci1);
                Yp = R(1,2)*X + R(2,2)*Y + gmfit.mu(i,ci2);
                
                plot(Xp, Yp,'LineWidth',2,'Color','k')
            end
            
            hold off
        end
    end
end