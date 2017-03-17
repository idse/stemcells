clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/')); 
dataDir = '/Volumes/IdseData/170210_ANBMP4withSB/part2';

%% metadata
%----------------

% positions are called Track and stored in separate directories
tracks = dir(fullfile(dataDir,'Track*'));
files = dir(fullfile(dataDir, tracks(1).name, '*oif'));
fname = fullfile(dataDir, tracks(1).name, files(1).name);

meta = Metadata(fname);
meta.nTime = numel(files);
meta.nPositions = numel(tracks);

% manual metadata
%-------------------

% TODO : modify MetadataAndor to contain all info below

treatmentTime = 0;
meta.timeInterval = '20 min';

meta.nWells = 4;
meta.posPerCondition = 4;
meta.conditions = {'B0','B1','B3','B10'};

% SET THIS TO TRUE IF MAKING AN '8-well' LOOP THROUGH A 4-WELL
loop4well = false;

nucChannel = 2;
S4Channel = 1;
tmax = meta.nTime;

%% save stitched previews of the MIPs

% TODO: remove black bands in between (make minimal value)
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
                    'tMax',             meta.nTime,...
                    'nucShrinkage',     1,...
                    'cytoSize',         6,...
                    'bgMargin',         10,...
                    'NCRcutoff',        [3 Inf]);

opts.cleanupOptions = struct('separateFused', true,...
    'clearBorder',true, 'minAreaStd', 1, 'minSolidity',0, 'minArea',300);


%% check that the options are set right

pi = 1;
P = Position(meta.nChannels, fname, meta.nTime);
P.setID(pi);
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

%%
tic

positions(meta.nPositions) = Position();

for pi = 1:meta.nPositions

    files = dir(fullfile(dataDir, tracks(pi).name, '*oif'));
    fname = fullfile(dataDir, tracks(pi).name, files(pi).name);

    positions(pi) = Position(meta.nChannels, fname, meta.nTime);
    positions(pi).setID(pi);
    positions(pi).extractData(dataDir, nucChannel, opts);
    positions(pi).makeTimeTraces();
    
    save(fullfile(dataDir,'positions'), 'positions');
end
toc

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

nWells = meta.nWells;
posPerCondition = meta.posPerCondition;

s = strsplit(meta.timeInterval,' ');
dt = str2double(s{1});
unit = s{2};
t = ((1:tmax) - treatmentTime)*dt;
axislim = [t(1), t(end)+50, 0.3 1.5];

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

for wellnr = 1%:nWells

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
    W = ones(size(W)); % don't weigh
    %W = bsxfun(@rdivide, posPerCondition*W, sum(W,2));
  
    nucTrace = cat(2,ttraceCat.nucLevelAvg);
    nucnucTrace = nucTrace(:,nucChannel:numel(opts.dataChannels):end);
    nucTrace = nucTrace(:,S4Channel:numel(opts.dataChannels):end);

    cytTrace = cat(2,ttraceCat.cytLevelAvg);
    cytTrace = cytTrace(:,S4Channel:numel(opts.dataChannels):end);
    
    bgTrace = cat(2,ttraceCat.background);
    bgTrace = bgTrace(:,S4Channel:numel(opts.dataChannels):end);
       
    ncellTrace = cat(1,positions.ncells)';
    ncellTrace = ncellTrace(:,conditionPositions);
    
    % determine what is bad to exclude it before averaging
    % not enough cells (complete focus failure or massive death):
    bad = ncellTrace < minNCells;
    
    % nuclear intensity differes too much from median
    % (temporary focal drift)
    medcut = 0.05;
    medNucTrace = medfilt2(nucnucTrace,[8 1]);
    badness = abs(medNucTrace - nucnucTrace)./abs(nucnucTrace);
    alsobad = badness > medcut;
    bad = bad | alsobad;
    
    nucTrace(bad) = NaN;
    cytTrace(bad) = NaN;
    bgTrace(bad) = NaN;
    ncellTrace(bad) = NaN;
    
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
            ratio(pi, bad(:,i)) = NaN;
            plot(t,ratio(pi,:),'Color', colors(i,:))
        end
        
        plot(t, ratioMean, graphfgc,'LineWidth',2)
        %plot(t, meanRatio, 'g','LineWidth',2)
        
        legend({'1','2','3','4','mean'},'Location','SouthEast');
        
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
            saveas(gcf, ['timeTrace_well' num2str(wellnr) '.png']);
        end
        
        hold off
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
yrange = [0.3, 1.4];

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
    W = ones(size(W)); % don't weigh
    %W = bsxfun(@rdivide, posPerCondition*W, sum(W,2));

    nucTrace = cat(2,ttraceCat.nucLevelAvg);
    nucTrace = nucTrace(:,S4Channel:numel(opts.dataChannels):end);
    cytTrace = cat(2,ttraceCat.cytLevelAvg);
    cytTrace = cytTrace(:,S4Channel:numel(opts.dataChannels):end);
    bgTrace = cat(2,ttraceCat.background);
    bgTrace = bgTrace(:,S4Channel:numel(opts.dataChannels):end);
    
    % determine what is bad to exclude it before averaging
    % not enough cells (complete focus failure or massive death):
    bad = ncellTrace < minNCells;
     
    % nuclear intensity differes too much from median
    % (temporary focal drift)
    medcut = 0.5;
    medNucTrace = medfilt2(nucnucTrace,[8 1]);
    badness = abs(medNucTrace - nucnucTrace)./abs(nucnucTrace);
    alsobad = badness > medcut;
    bad = bad | alsobad;
    
    nucTrace(bad) = NaN;
    cytTrace(bad) = NaN;
    bgTrace(bad) = NaN;
    ncellTrace(bad) = NaN;
    
    nucMean = nanmean(nucTrace(1:tmax,:).*W,2);
    cytMean = nanmean(cytTrace(1:tmax,:).*W,2);
    bgMean = nanmean(bgTrace(1:tmax,:).*W,2);
    
    ratioMean = (nucMean-bgMean)./(cytMean - bgMean);
    
    %bad = any(cat(1,positions(conditionPositions).ncells) < minNCells,1);
    %ratioMean(any(bad,2)) = NaN;

    %ratioMean = ratioMean -  baseline(wellnr) + mean(baseline);
    plot(t, ratioMean,'LineWidth',2,'Color',colors(wellidx,:))
    ylim(yrange);

    %plot(t, bgMean,'LineWidth',2,'Color',colors(wellidx,:))
    %plot(t, nucMean,'LineWidth',2,'Color',colors(wellidx,:))
    %plot(t, cytMean,'LineWidth',2,'Color',colors(wellidx,:))
    xlim([t(1), t(end)+50]);

    fs = 24;
    xlabel(['time (' unit ')'], 'FontSize',fs,'FontWeight','Bold')
    ylabel('nuclear : cytoplasmic Smad4', 'FontSize',fs,'FontWeight','Bold');

    set(gcf,'color','w');
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
end
hold off

set(gca,'FontSize', 16)
legend(meta.conditions(wellsWanted),'Location','SouthEast');
if saveResult
    %export_fig(['timeTrace_multipleConditions.pdf'],'-native -m2');
    saveas(gcf,'timeTrace_multipleConditions.fig');
    saveas(gcf,'timeTrace_multipleConditions.png');
end
