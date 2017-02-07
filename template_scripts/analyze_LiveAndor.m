clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
dataDir = '/Users/idse/data_tmp/170120_cyclohexRep';

meta = MetadataAndor(dataDir);

% manual metadata
%-------------------

% TODO : modify MetadataAndor to contain all info below

treatmentTime = 4;

meta.nWells = 6;
meta.posPerCondition = 4;
meta.conditions = {'SB+M','A+M+C','A','Nog+M','A+M','A+C-1hr-SB','A+C','M'};

% SET THIS TO TRUE IF MAKING AN '8-well' LOOP THROUGH A 4-WELL
loop4well = true;

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

opts = struct(...
                    'dataChannels',     [S4Channel nucChannel],...
                    'fgChannel',        S4Channel,...
                    'cytoplasmicLevels',true,... 
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'MIPidxDir',        fullfile(dataDir,'MIP'),...
                    'tMax',             5,...
                    'nucShrinkage',     2,...
                    'cytoSize',         8,...
                    'bgMargin',         10,...
                    'NCRcutoff',        [3 Inf]);

opts.cleanupOptions = struct('separateFused', true,...
    'clearBorder',true, 'minAreaStd', 1, 'minSolidity',0, 'minArea',500);


%% check that the options are set right

pi = 1;
P = DynamicPositionAndor(meta, pi);
time = 1;
opts.tMax = time;

% try out the nuclear cleanup settings on some frame:
% bla = nuclearCleanup(seg(:,:,time), opts.cleanupOptions);
% imshow(bla)

debugInfo = P.extractData(dataDir, nucChannel, opts);

bgmask = debugInfo.bgmask;
nucmask = debugInfo.nucmask;
cytmask = false(size(nucmask));
cytmask(cat(1,debugInfo.cytCC.PixelIdxList{:}))=true;

bg = P.cellData(time).background
nucl = P.cellData(time).nucLevelAvg
cytl = P.cellData(time).cytLevelAvg
(nucl-bg)./(cytl - bg)

im = P.loadImage(dataDir, S4Channel, time);
MIP = max(im,[],3);
A = imadjust(mat2gray(MIP));
s = 0.4;
imshow(cat(3, A + 0*bgmask, A + s*nucmask, A + s*cytmask));

opts.tMax = tmax;

%% run the analysis on all time points

tic
positions(meta.nPositions) = DynamicPositionAndor();

for pi = 1:meta.nPositions

    positions(pi) = DynamicPositionAndor(meta, pi);
    positions(pi).extractData(dataDir, nucChannel, opts);
    positions(pi).makeTimeTraces();
    save(fullfile(dataDir,'positions'), 'positions');
end
toc

%% load results if above block was run previously

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

nWells = meta.nWells;
posPerCondition = meta.posPerCondition;

s = strsplit(meta.timeInterval,' ');
dt = str2double(s{1});
unit = s{2};
t = ((1:tmax) - treatmentTime)*dt;
axislim = [t(1), t(end)+50, 0.3, 1.3];

frame = {};
cd(dataDir);
saveResult = true; % CHECK
minNCells = 10; % minimal number of cells
fgc = 'k';
bgc = 'w';
graphbgc = 1*[1 1 1]; 
graphfgc = 'r';
%w, k, 0.5, w

colors = hsv(posPerCondition); %0.5*[1 0 0]

baseline = zeros([1 nWells]); % store baseline avg of each well

for wellnr = 4%:nWells

    if loop4well
        flipwell = 9-wellnr;
        conditionPositions = [(posPerCondition/2)*(wellnr-1)+1:(posPerCondition/2)*wellnr...
                              (posPerCondition/2)*(flipwell-1)+1:(posPerCondition/2)*flipwell];
    else
        conditionPositions = posPerCondition*(wellnr-1)+1:posPerCondition*wellnr;
    end
    
    ttraceCat = cat(1,positions.timeTraces);
    ttraceCat = ttraceCat(conditionPositions);

    % weigh by number of cells, tends to make little difference
    W = cat(1,positions.ncells); 
    W = W(conditionPositions,:)';
    %W = ones(size(W)); % don't weigh
    W = bsxfun(@rdivide, posPerCondition*W, sum(W,2));
    
    nucTrace = cat(2,ttraceCat.nucLevelAvg);
    cytTrace = cat(2,ttraceCat.cytLevelAvg);
    bgTrace = cat(2,ttraceCat.background);
    
    nucMean = nanmean(nucTrace(1:tmax,:).*W,2);
    cytMean = nanmean(cytTrace(1:tmax,:).*W,2);
    bgMean = nanmean(bgTrace(1:tmax,:).*W,2);
    
    ratioMean = (nucMean - bgMean)./(cytMean - bgMean);
    %meanRatio = nanmean(ratio,1); % makes little difference
    baseline(wellnr) = mean(ratioMean(t < treatmentTime));
    
    % THIS SHOULD BE REARRANGED WITH ti ON THE INSIDE AND pi OUTSIDE
    
    for ti = 1%:positions(1).nTime
        clf
        hold on
        ratio = zeros([numel(conditionPositions) tmax]);
        for i = 1:numel(conditionPositions)

            pi = conditionPositions(i);
            
            nucTrace = positions(pi).timeTraces.nucLevelAvg;
            bgTrace = positions(pi).timeTraces.background;
            cytTrace = positions(pi).timeTraces.cytLevelAvg;
            
            R = (nucTrace - bgTrace)./(cytTrace - bgTrace);
            ratio(pi,:) = R(1:tmax)';
            ratio(pi, positions(pi).ncells < minNCells) = NaN;
            plot(t,ratio(pi,:),'Color', colors(i,:))
        end
        
        plot(t, ratioMean, graphfgc,'LineWidth',2)
        %plot(t, meanRatio, 'g','LineWidth',2)
        
        %legend({'1','2','3','4','mean'});
        
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
            export_fig(['timeTrace_well' num2str(wellnr) '.png'],'-native -m2');
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
figure,
s = strsplit(meta.timeInterval,' ');
dt = str2double(s{1});
unit = s{2};
t = ((1:tmax) - treatmentTime)*dt;

posPerCondition = meta.posPerCondition;

frame = {};
cd(dataDir);
saveResult = true;
minNCells = 10; % minimal number of cells

wellsWanted = 1:4;
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
    if loop4well
        flipwell = 9-wellnr;
        conditionPositions = [(posPerCondition/2)*(wellnr-1)+1:(posPerCondition/2)*wellnr...
                              (posPerCondition/2)*(flipwell-1)+1:(posPerCondition/2)*flipwell];
    else
        conditionPositions = posPerCondition*(wellnr-1)+1:posPerCondition*wellnr;
    end
    
    ttraceCat = cat(1,positions.timeTraces);
    ttraceCat = ttraceCat(conditionPositions);
    
    % weigh by number of cells
    W = cat(1,positions.ncells); 
    W = W(conditionPositions,:)';
    %W = ones(size(W)); % don't weigh
    W = bsxfun(@rdivide, posPerCondition*W, sum(W,2));

    nucTrace = cat(2,ttraceCat.nucLevelAvg);
    cytTrace = cat(2,ttraceCat.cytLevelAvg);
    bgTrace = cat(2,ttraceCat.background);
    
    nucMean = nanmean(nucTrace(1:tmax,:).*W,2);
    cytMean = nanmean(cytTrace(1:tmax,:).*W,2);
    bgMean = nanmean(bgTrace(1:tmax,:).*W,2);
    
    ratioMean = (nucMean - bgMean)./(cytMean - bgMean);
    
%     bad = any(cat(1,positions(conditionPositions).ncells) < minNCells,1);
%     ratioMean(bad) = NaN;

    %ratioMean = ratioMean -  baseline(wellnr) + mean(baseline);
    plot(t, ratioMean,'LineWidth',2,'Color',colors(wellidx,:))
    ylim([0.3, 1.8]);

    %plot(t, bgMean,'LineWidth',2,'Color',colors(wellidx,:))
    %plot(t, nucMean,'LineWidth',2,'Color',colors(wellidx,:))
    %plot(t, cytMean,'LineWidth',2,'Color',colors(wellidx,:))
    xlim([t(1), t(end)+50]);

    g0 = [0 1 0];
    if max(ratioMean) > 1
        g0 = g0/max(ratioMean);
    end

    fs = 24;
    xlabel(['time (' unit ')'], 'FontSize',fs,'FontWeight','Bold')
    ylabel('nuclear : cytoplasmic Smad4', 'FontSize',fs,'FontWeight','Bold');

    set(gcf,'color','w');
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
end
hold off
title('washing protocols for A10 pulses');
set(gca,'FontSize', 16)
legend(meta.conditions(wellsWanted));
if saveResult
    export_fig(['timeTrace_multipleConditions.pdf'],'-native -m2');
    saveas(gcf,'timeTrace_multipleConditions.fig');
end

%% combine ratio statistics of positions in each well

pi = 1;
ti = 1;
wellsWanted = 1:meta.nWells;
cellDataWells = cell([numel(wellsWanted) meta.nTime]);

for wellidx = 1:numel(wellsWanted)
    
    wellnr = wellsWanted(wellidx);
    conditionPositions = meta.posPerCondition*(wellnr-1)+1:meta.posPerCondition*wellnr;
    
    for i = 1:numel(conditionPositions)

        for ti = 1:meta.nTime
            pi = conditionPositions(i);    
            CD = positions(pi).cellData(ti);
            ratio = (CD.nucLevel - CD.background)./(CD.cytLevel - CD.background); 
            cellDataWells{wellidx, ti} = cat(1,cellDataWells{wellidx,ti}, ratio);
        end
    end
end

%% time evolution of distribution in some well

wi = 1;

rmin = 0;
rmax = 3;
bins = linspace(rmin,rmax,50);
clf 
tidx = 1:4:80;
colors = hsv(numel(tidx));
hold on
for ti = 1:numel(tidx)
    n=histc(cellDataWells{wi,tidx(ti)}, bins);
    plot(bins,n,'Color',colors(ti,:));
end
hold off
xlim([rmin rmax]);

%% compare late time distribution in different wells

figure,
ti = 1;

rmin = 0;
rmax = 3;
bins = linspace(rmin,rmax,50);

colors = [0 0 1; 0 0.5 0; 1 0 0; 1 0 0; 0 0.5 0; 0 0 1];
clf 
hold on

for wi = 1:meta.nWells

    n=histc(cellDataWells{wi,ti}, bins);
    n = n./sum(n);
    n = cumsum(n);
    plot(bins,n,'Color', colors(wi,:));
end
hold off
xlim([rmin rmax]);
ylim([0 1]);

legend(meta.conditions);

%% integrals

treatmentTime = 10;
t = ((1:tmax) - treatmentTime)*dt;

lateT = treatmentTime + 420/dt;

F = ratioMean - min(ratioMean);
I = cumsum(F*dt);

Ipulse = I;
Ipulse(lateT:end) = Ipulse(lateT);
Ibase = mean(F(lateT:end))*t;

thr = t/60;
plot(thr,I)
hold on
plot(thr,Ibase,'r')
plot(thr,Ipulse,'g')
hold off
xlim([thr(1) thr(end)]);
xlabel('time (hours)');

legend({'adaptive response','pulse','baseline'},'Location','SouthEast');
title('integrated signal');
