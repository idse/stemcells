clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
dataDir = '/Users/idse/data_tmp/160812_siRNASki+Skil';

meta = MetadataAndor(dataDir);

% manual metadata
%-------------------

% TODO : modify MetadataAndor to contain all info below

%barefname = '160812_siRNASki+Skil';
barefname = '160812';
treatmentTime = 4;

meta.nWells = 6;
meta.posPerCondition = 4;
meta.conditions = {'ctrl','siSmad4', 'siSkiSnoN', 'siSkiSnoN',...
              'siSmad4','ctrl'};

% barefname = 'SBbackground';
% treatmentTime = 4;
% posPerCondition = 4;
% nWells = 8;

nucChannel = 2;
S4Channel = 1;
tmax = meta.nTime;

% visualize positions
%---------------------

% meta.displayPositions;

% TODO: create merged cellData for montage
% movies of distribution over time

           
%% save stitched previews of the MIPs

stitchedPreviews(dataDir, meta); 

%% extract nuclear and cytoplasmic levels

% externally I will have all indices starting at 1
% the Andor offset to start at 0 will be internal

opts = struct(  'cytoplasmicLevels',    true,... %'tMax', 25,...
                    'dataChannels',     S4Channel,...
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'nucShrinkage',     2,...
                    'MIPidxDir',        fullfile(dataDir,'MIP'),...
                    'tMax',             tmax);

opts.cleanupOptions = struct('separateFused', true,...
    'clearBorder',true, 'minAreaStd', 1, 'minSolidity',0, 'minArea',1500);

% try out the setting on some frame:
% bla = nuclearCleanup(seg(:,:,50), opts.cleanupOptions);
% imshow(bla)

tic
positions(meta.nPositions) = DynamicPositionAndor();

for pi = 1:meta.nPositions

    positions(pi) = DynamicPositionAndor(meta, pi);
    positions(pi).extractData(dataDir, nucChannel, opts);
    positions(pi).makeTimeTraces();
end
toc

%save(fullfile(dataDir,'positions'), 'positions');
%['positions_' datestr(now,'yymmdd')]

%%
load(fullfile(dataDir,'positions'));

% positions_before = load(fullfile(dataDir,'positions_before'));
% positions_before = positions_before.positions;
% 
% for pi = 1:meta.nPositions
%     positions(pi).cellData = cat(2, positions_before(pi).cellData,...
%                                             positions(pi).cellData);
%     positions(pi).ncells = cat(2, positions_before(pi).ncells,...
%                                             positions(pi).ncells);
% 	positions(pi).nTime = positions_before(pi).nTime + positions(pi).nTime;
%     positions(pi).makeTimeTraces;
% end
% 
% save(fullfile(dataDir,'positions'), 'positions');

%% make a video (and figure) of the time traces

s = strsplit(meta.timeInterval,' ');
dt = str2double(s{1});
unit = s{2};
t = ((1:tmax) - treatmentTime)*dt;
axislim = [t(1), t(end)+50, 0.4, 1.6];

frame = {};
cd(dataDir);
saveResult = true; % CHECK
minNCells = 10; % minimal number of cells
fgc = 'k';
bgc = 'w';
graphbgc = 1*[1 1 1]; 
graphfgc = 'r';
%w, k, 0.5, w

baseline = zeros([1 nWells]); % store baseline avg of each well

for wellnr = 1:nWells

    conditionPositions = posPerCondition*(wellnr-1)+1:posPerCondition*wellnr;

    ttraceCat = cat(1,positions.timeTraces);
    ttraceCat = ttraceCat(conditionPositions);
    
    % weigh by number of cells, tends to make little difference
    W = cat(1,positions.ncells); 
    W = W(conditionPositions,:)';
    %W = ones(size(W)); % don't weigh
    
    nucTrace = cat(2,ttraceCat.nucLevelAvg);
    cytTrace = cat(2,ttraceCat.cytLevelAvg);
    bgTrace = cat(2,ttraceCat.background);
    
    nucMean = nanmean(nucTrace(1:tmax,:).*W,2)./sum(W,2);
    cytMean = nanmean(cytTrace(1:tmax,:).*W,2)./sum(W,2);
    bgMean = nanmean(bgTrace(1:tmax,:).*W,2)./sum(W,2);
    
    ratioMean = (nucMean - bgMean)./(cytMean - bgMean);
    %meanRatio = nanmean(ratio,1); % makes little difference
    baseline(wellnr) = mean(ratioMean(t < treatmentTime));
        
    % THIS SHOULD BE REARRANGED WITH ti ON THE INSIDE AND pi OUTSIDE
    
    for ti = 1%:positions(1).nTime
        clf 
        hold on
        ratio = zeros([numel(conditionPositions) tmax]);
        for pi = conditionPositions

            nucTrace = positions(pi).timeTraces.nucLevelAvg;
            bgTrace = positions(pi).timeTraces.background;
            cytTrace = positions(pi).timeTraces.cytLevelAvg;
            
            R = (nucTrace - bgTrace)./(cytTrace - bgTrace);
            ratio(pi,:) = R(1:tmax)';
            ratio(pi, positions(pi).ncells < minNCells) = NaN;
            plot(t,ratio(pi,:),'Color', 0.5*[1 0 0])
        end
        
        plot(t, ratioMean, graphfgc,'LineWidth',2)
        %plot(t, meanRatio, 'g','LineWidth',2)
        
        fs = 24;
        xlabel(['time (' unit ')'], 'FontSize',fs,'FontWeight','Bold','Color',fgc)
        ylabel('nuclear : cytoplasmic Smad4', 'FontSize',fs,'FontWeight','Bold','Color',fgc);
        
        axis(axislim);
        set(gcf,'color',bgc);
        set(gca, 'LineWidth', 2);
        set(gca,'FontSize', fs)
        set(gca,'FontWeight', 'bold')
        set(gca,'XColor',fgc);
        set(gca,'YColor',fgc);
        set(gca,'Color',graphbgc);
        
        if saveResult
            export_fig(['timeTrace_well' num2str(wellnr) '.pdf'],'-native -m2');
        end

        if ~isnan(ratioMean(ti))
            g0 = [0 1 0];
            if max(ratioMean) > 1
                g0 = g0/max(ratioMean);
            end
            plot(t(ti), ratioMean(ti),'o','LineWidth',3,'MarkerEdgeColor',g0,...
                                    'MarkerSize',20,'MarkerFaceColor',ratioMean(ti)*g0)
        end
        hold off

        frame{ti} = export_fig(gcf,'-native -m2');
    end
%     if saveResult
%         disp('saving');
%         v = VideoWriter(fullfile(dataDir,['ratioplot_white_well' num2str(wellnr) '.mp4']),'MPEG-4');
%         v.FrameRate = 5;
%         open(v)
%         for ti = 1:positions(1).nTime
%             writeVideo(v,frame{ti})
%         end
%         close(v);
%     end
end

%% make a combined plot of several conditions

s = strsplit(meta.timeInterval,' ');
dt = str2double(s{1});
unit = s{2};
t = ((1:tmax) - treatmentTime)*dt;

frame = {};
cd(dataDir);
saveResult = true;
minNCells = 10; % minimal number of cells

wellsWanted = [1 2 3];
%wellsWanted = 4:8;
colors = lines(numel(wellsWanted));

clf 
hold on

for wellidx = 1:numel(wellsWanted)
    
    wellnr = wellsWanted(wellidx);

%     if wellnr == 2
%         conditionPositions = 4*(wellnr-1)+1;
%     else
%         conditionPositions = 4*(wellnr-1)+1:4*wellnr;
%     end

    conditionPositions = posPerCondition*(wellnr-1)+1:posPerCondition*wellnr;
    
    ttraceCat = cat(1,positions.timeTraces);
    ttraceCat = ttraceCat(conditionPositions);
    
    % weigh by number of cells
    W = cat(1,positions.ncells); 
    W = W(conditionPositions,:)';
    %W = ones(size(W)); % don't weigh
    
    nucTrace = cat(2,ttraceCat.nucLevelAvg);
    cytTrace = cat(2,ttraceCat.cytLevelAvg);
    bgTrace = cat(2,ttraceCat.background);
    
    nucMean = nanmean(nucTrace(1:tmax,:).*W,2)./sum(W,2);
    cytMean = nanmean(cytTrace(1:tmax,:).*W,2)./sum(W,2);
    bgMean = nanmean(bgTrace(1:tmax,:).*W,2)./sum(W,2);
    
    ratioMean = (nucMean - bgMean)./(cytMean - bgMean);
    
%     bad = any(cat(1,positions(conditionPositions).ncells) < minNCells,1);
%     ratioMean(bad) = NaN;

    %ratioMean = ratioMean -  baseline(wellnr) + mean(baseline);
    plot(t, ratioMean,'LineWidth',2,'Color',colors(wellidx,:))

    g0 = [0 1 0];
    if max(ratioMean) > 1
        g0 = g0/max(ratioMean);
    end

    fs = 24;
    xlabel(['time (' unit ')'], 'FontSize',fs,'FontWeight','Bold')
    ylabel('nuclear : cytoplasmic Smad4', 'FontSize',fs,'FontWeight','Bold');

    axis([t(1), t(end)+50, 0.6, 1.6]);
    set(gcf,'color','w');
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
end
hold off
%title('comparison');
set(gca,'FontSize', 16)
legend(conditions(wellsWanted));
if saveResult
    export_fig(['timeTrace_multipleConditions2.pdf'],'-native -m2');
end