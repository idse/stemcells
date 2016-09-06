clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/')); 
dataDir = '/Volumes/IdseData3/160902_C2C12pluronic';

list = dir(fullfile(dataDir,'*vsi'));
vsifile = fullfile(dataDir,list(1).name);

%% metadata
%----------------

metaDataFile = fullfile(dataDir,'metaData.mat');

% read metadata
if exist(metaDataFile,'file')
    disp('loading previously stored metadata');
    load(metaDataFile);
elseif exist('vsifile','var') && exist(vsifile,'file')
    disp('extracting metadata from vsi file');
    meta = Metadata(vsifile);
end

% manually entered metadata
%------------------------------

meta.timeInterval = '9 min';
meta.nPositions = 16;

%save(fullfile(dataDir,'metaData'),'meta');

barefname = 'syringeC2C12_7hintervals';
treatmentTime = 7;
meta.nWells = 2;
meta.posPerCondition = 6;
meta.conditions = {'pluronic + fibronectin', 'fibronectin'};

nucChannel = 2;
S4Channel = 1;
           
%% save stitched previews of the MIPs

% TODO: remove black bands in between (make minimal value)
stitchedPreviews(dataDir, meta);

%% extract nuclear and cytoplasmic levels

% externally I will have all indices starting at 1
% the Andor offset to start at 0 will be internal

opts = struct(  'cytoplasmicLevels',    true,... %'tMax', 25,...
                    'dataChannels',     S4Channel,...
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'nucShrinkage',     2);
                    
% add this if there is a MIPidx
                    %'MIPidxDir',        fullfile(dataDir,'MIP'));

opts.cleanupOptions = struct('separateFused', true,...
    'clearBorder',true, 'minAreaStd', 1, 'minSolidity',0, 'minArea',800);

%%
% try out the setting on some frame:
pos = DynamicPositionAndor(meta,1);
seg = pos.loadSegmentation(fullfile(dataDir,'MIP'),nucChannel);
bla = nuclearCleanup(seg(:,:,1), opts.cleanupOptions);
imshow(bla)

%%
tic
%positions(meta.nPositions) = DynamicPositionAndor();
for pi = 1:meta.nPositions

    positions(pi) = DynamicPositionAndor(meta,pi);
    positions(pi).extractData(fullfile(dataDir,'MIP'), nucChannel, opts);
    positions(pi).makeTimeTraces();
end
toc

save(fullfile(dataDir,'positions'), 'positions');
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

%% make a video of the time traces

s = strsplit(meta.timeInterval,' ');
dt = str2double(s{1});
unit = s{2};
t = ((1:positions(1).nTime) - treatmentTime)*dt;

frame = {};
cd(dataDir);
saveResult = true;
minNCells = 10; % minimal number of cells
fgc = 'w';
bgc = 'k';

for wellnr = 1:nWells

    conditionPositions = posPerCondition*(wellnr-1)+1:posPerCondition*wellnr;

    ttraceCat = cat(1,positions.timeTraces);
    ttraceCat = ttraceCat(conditionPositions);
    nucMean = mean(cat(2,ttraceCat.nucLevelAvg),2);
    cytMean = mean(cat(2,ttraceCat.cytLevelAvg),2);
    bgMean = mean(cat(2,ttraceCat.background),2);

    % THIS SHOULD BE REARRANGED WITH ti ON THE INSIDE AND pi OUTSIDE
    
    for ti = 1:positions(1).nTime
        clf 
        hold on
        for pi = conditionPositions

            nucTrace = positions(pi).timeTraces.nucLevelAvg;
            bgTrace = positions(pi).timeTraces.background;
            cytTrace = positions(pi).timeTraces.cytLevelAvg;

            ratio = (nucTrace - bgTrace)./(cytTrace - bgTrace);
            ratio(positions(pi).ncells < minNCells) = NaN;
            plot(t,ratio,'Color', 0.5*[1 0 0])
        end
        ratioMean = (nucMean - bgMean)./(cytMean - bgMean);
        bad = any(cat(1,positions(conditionPositions).ncells) < minNCells,1);
        ratioMean(bad) = NaN;
        plot(t, ratioMean,'w','LineWidth',2)
        
        fs = 24;
        xlabel(['time (' unit ')'], 'FontSize',fs,'FontWeight','Bold','Color',fgc)
        ylabel('nuclear : cytoplasmic Smad4', 'FontSize',fs,'FontWeight','Bold','Color',fgc);
        
        axis([t(1), t(end)+50, 0.3, 2]);
        set(gcf,'color',bgc);
        set(gca, 'LineWidth', 2);
        set(gca,'FontSize', fs)
        set(gca,'FontWeight', 'bold')
        set(gca,'XColor',fgc);
        set(gca,'YColor',fgc);
        set(gca,'Color',0.5*[1 1 1]);
        
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
    if saveResult
        v = VideoWriter(fullfile(dataDir,['ratioplot_well' num2str(wellnr) '.mp4']),'MPEG-4');
        v.FrameRate = 5;
        open(v)
        for ti = 1:positions(1).nTime
            writeVideo(v,frame{ti})
        end
        close(v);
    end
end

%% make a combined plot of several conditions

s = strsplit(meta.timeInterval,' ');
dt = str2double(s{1});
unit = s{2};
t = ((1:positions(1).nTime) - treatmentTime)*dt;

frame = {};
cd(dataDir);
saveResult = true;
minNCells = 10; % minimal number of cells

wellsWanted = 1;
colors = lines(numel(wellsWanted));

% pulse graph information
plotPulseGraph = true;
logDir = dataDir;
logName = '160615_144240_timeLog.txt';
startTime = treatmentTime*dt;

clf 
hold on

for wellidx = 1:numel(wellsWanted)
    
    wellnr = wellsWanted(wellidx);

    % find the positions for wellnr
    conditionPositions = posPerCondition*(wellnr-1)+1:posPerCondition*wellnr;

    % retrieve fluorescent data
    ttraceCat = cat(1,positions.timeTraces);
    ttraceCat = ttraceCat(conditionPositions);
    
    % these arrays have zeros where NaNs would be
    nuc = cat(2,ttraceCat.nucLevelAvg);
    cyt = cat(2,ttraceCat.cytLevelAvg);
    bg = cat(2,ttraceCat.background);
    
    % the zeros in these arrays are changed to NaN
    if wellidx == 1
        warning('All fluorescent levels of zero are assumed to be misimaged.');
    end
    nuc(nuc == 0) = NaN;
    cyt(cyt == 0) = NaN; 
    bg(bg == 0) = NaN;
    
    % taking the ratios for each position
    ratios = zeros(size(nuc));
    for pi = 1:posPerCondition
        ratios(:,pi) = (nuc(:,pi) - bg(:,pi))./(cyt(:,pi) - bg(:,pi));
        % for this position, make time points with too few cells = NaN
        bad = cat(1,positions(posPerCondition*(wellnr-1) + pi).ncells) < minNCells;
        ratios(bad',pi) = NaN;
    end
    
    % plot the average of the ratios
    ratioMean = mean(ratios,2,'omitnan');
    plot(t, ratioMean,'LineWidth',2,'Color',colors(wellidx,:))

    g0 = [0 1 0];
    if max(ratioMean) > -7
        g0 = g0/max(ratioMean);
    end
    
    % plot pulse graph
    if wellidx == 1
        tmp = ratioMean;
    else
        tmp = [tmp ratioMean];
    end
    if plotPulseGraph && wellidx == numel(wellsWanted)
        plotPulse(logDir,logName,startTime,min(tmp(:)),max(tmp(:)))
    end

    % plot parameters
    fs = 24;
    xlabel(['time (' unit ')'], 'FontSize',fs,'FontWeight','Bold')
    ylabel('nuclear : cytoplasmic Smad4', 'FontSize',fs,'FontWeight','Bold');

    axis([t(1), t(end)+50, 0.3, 2]);
    set(gcf,'color','w');
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')

    %frame{ti} = export_fig(gcf,'-native -m2');
end
hold off
title('Syringe Pump: C2C12 w/ tgfb @ 7h intervals');
set(gca,'FontSize', 16)
legend(conditions(wellsWanted));
if saveResult
    export_fig(['timeTrace_multipleConditions.pdf'],'-native -m2');
end