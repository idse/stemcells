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

% determine process numbers of vsi files
vsifiles = dir(fullfile(dataDir,'*.vsi'));
processNumbers = zeros([numel(vsifiles) 1],'uint16');
for i = 1:numel(vsifiles)
    s = strsplit(vsifiles(i).name,{'_','.'});
    processNumbers(i) = uint16(str2double(s{2}));
end
processnr = min(processNumbers):max(processNumbers);

% manually entered metadata
%------------------------------

meta.timeInterval = '9 min';

%save(fullfile(dataDir,'metaData'),'meta');

barefname = 'syringeC2C12_7hintervals';
treatmentTime = 7;
meta.nWells = 2;
meta.posPerCondition = 6;
meta.nPositions = meta.nWells*meta.posPerCondition;
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
                    'fgChannel',        S4Channel,...
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'tMax',             meta.nTime,...
                    'nucShrinkage',     2,...
                    'cytoSize',         10,...
                    'bgMargin',         4);
                
% add this if there is a MIPidx
                    %'MIPidxDir',        fullfile(dataDir,'MIP'));

opts.cleanupOptions = struct('separateFused', true,...
    'clearBorder',true, 'minAreaStd', 1, 'minSolidity',0, 'minArea',800);

%%
tic
positions(meta.nPositions) = Position();
for pi = 1:meta.nPositions

    vsifile = fullfile(dataDir,['Process_' num2str(processnr(pi)) '.vsi']);

    positions(pi) = Position(meta.nChannels, vsifile, meta.nTime);
    positions(pi).setID(pi);
    positions(pi).extractData(dataDir, nucChannel, opts);
    positions(pi).makeTimeTraces();
    
    save(fullfile(dataDir,'positions'), 'positions');
end
toc

%% finding the right options for nuclear cleanup

pi = 1;
vsifile = fullfile(dataDir,['Process_' num2str(processnr(pi)) '.vsi']);
P = Position(meta.nChannels, vsifile, meta.nTime);
P.setID(pi);
seg = P.loadSegmentation(fullfile(dataDir,'MIP'), nucChannel);

% try out the nuclear cleanup settings on some frame:
time = 1;
bla = nuclearCleanup(seg(:,:,time), opts.cleanupOptions);
figure, imshow(cat(3,mat2gray(bla),seg(:,:,time),seg(:,:,time)))

%% finding the right options for extract data

opts.tMax = time;
debugInfo = P.extractData(dataDir, nucChannel, opts);

bgmask = debugInfo.bgmask;
nucmask = debugInfo.nucmask;
cytmask = false(size(nucmask));
cytmask(cat(1,debugInfo.cytCC.PixelIdxList{:}))=true;

bg = P.cellData(time).background
nucl = P.cellData(time).nucLevelAvg
cytl = P.cellData(time).cytLevelAvg
(nucl-bg)/(cytl - bg)

im = P.loadImage(dataDir, S4Channel, time);
MIP = max(im,[],3);
A = imadjust(mat2gray(MIP));
s = 0.4;
imshow(cat(3, A + 0*bgmask, A + s*nucmask, A + s*cytmask));

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

tmax = meta.nTime;
nWells = meta.nWells;
posPerCondition = meta.posPerCondition;

s = strsplit(meta.timeInterval,' ');
dt = str2double(s{1});
unit = s{2};
t = ((1:tmax) - treatmentTime)*dt;
axislim = [t(1), t(end)+50, 0.5, 1.6];

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

for wellnr = 1:nWells

    conditionPositions = posPerCondition*(wellnr-1)+1:posPerCondition*wellnr;

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
    figure
    for ti = 1%:positions(1).nTime
        
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
        
        legendstr = {};
        for i = 1:posPerCondition
            legendstr = [legendstr num2str(i)];
        end
        legendstr = [legendstr 'mean'];
        legend(legendstr);
        
        title(meta.conditions{wellnr});
        
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

%         if ~isnan(ratioMean(ti))
%             g0 = [0 1 0];
%             if max(ratioMean) > 1
%                 g0 = g0/max(ratioMean);
%             end
%             plot(t(ti), ratioMean(ti),'o','LineWidth',3,'MarkerEdgeColor',g0,...
%                                     'MarkerSize',20,'MarkerFaceColor',ratioMean(ti)*g0)
%         end
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

posPerCondition = meta.posPerCondition;

s = strsplit(meta.timeInterval,' ');
dt = str2double(s{1});
unit = s{2};
t = ((1:positions(1).nTime) - treatmentTime)*dt;

frame = {};
cd(dataDir);
saveResult = true;
minNCells = 10; % minimal number of cells

wellsWanted = 1:2;
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
    ratioMean = nanmean(ratios,2);
    plot(t, ratioMean,'LineWidth',2,'Color',colors(wellidx,:))

    g0 = [0 1 0];
    if max(ratioMean) > -7
        g0 = g0/max(ratioMean);
    end
    
%     % plot pulse graph
%     if wellidx == 1
%         tmp = ratioMean;
%     else
%         tmp = [tmp ratioMean];
%     end
%     if plotPulseGraph && wellidx == numel(wellsWanted)
%         plotPulse(logDir,logName,startTime,min(tmp(:)),max(tmp(:)))
%     end

    % plot parameters
    fs = 24;
    xlabel(['time (' unit ')'], 'FontSize',fs,'FontWeight','Bold')
    ylabel('nuclear : cytoplasmic Smad4', 'FontSize',fs,'FontWeight','Bold');

    axis([t(1), t(end)+50, 0.6, 1.5]);
    set(gcf,'color','w');
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')

    %frame{ti} = export_fig(gcf,'-native -m2');
end
hold off
%title('Syringe Pump: C2C12 w/ tgfb @ 7h intervals');
set(gca,'FontSize', 16)
legend(meta.conditions(wellsWanted));
if saveResult
    export_fig(['timeTrace_multipleConditions.pdf'],'-native -m2');
end