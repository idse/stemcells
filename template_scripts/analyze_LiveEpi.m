clear all; close all;

addpath(genpath('C:\Users\Thomas\Documents\GitHub\stemcells')); 

%dataDir = '/Volumes/Seagate Backup Plus Drive/160317_ibidi_RIactivin';
%dataDir = '/Users/idse/data_tmp/cycloheximide_after_20160330_32055 PM';
%dataDir = '/Users/idse/data_tmp/cycloheximide_before_20160330_42945 PM';
%dataDir = '/Volumes/IdseData/160416_RIvsnoRI';
dataDir = 'D:\160617_syringeC2C12_7hintervals';


vsifile = fullfile(dataDir,'Process_12.vsi');
maxMemoryGB = 4;


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
    
% or pretty much enter it by hand
else
    meta = Metadata();
    h = Tiff(btfname);
    meta.ySize = h.getTag('ImageLength');
    meta.xSize = h.getTag('ImageWidth');
    h.close;
    meta.nChannels = 4;
    meta.channelNames = {'DAPI','GFP','RFP','CY5'};
    meta.xres = 0.650/2;
    meta.yres = meta.xres;
end

% manually entered metadata
%------------------------------

meta.channelLabel = {'H2B','Smad4'};
%meta.channelLabel = {'Bra','H2B', 'Smad4', 'Sox17'};
%meta.tPerFile = 138;
meta.filename = 'syringeC2C12_7hintervals_MIP_p%.4d_w%.4d.tif';
meta.nTime = 200;
meta.timeInterval = '10 min';
meta.nPositions = 8;

%save(fullfile(dataDir,'metaData'),'meta');


% 
% meta = MetadataAndor(dataDir);
% %meta.nTime = 110; % JUST FOR cycloheximide_after
% %meta.nTime = 3; % JUST FOR cycloheximide_before
% filenameFormat = meta.filename;
% 
% % manual metadata
% %-------------------
% 
% % TODO : modify MetadataAndor to contain all info below
% 
% % barefname = 'RIactivin100';
% % treatmentTime = 8; % first time point after treatment
% % conditions = {'RI + Activin 100 ng/ml'};
% % posPerCondition = 16;
% % nWells = 1;
% 
% %barefname = 'cycloheximide_after';
% %barefname = 'cycloheximide_before';
% % treatmentTime = 4;
% % conditions = {'no treatment','Activin 100 ng/ml', 'BMP 50 ng/ml', 'cyclohex 50 \mu g/ml',...
% %               'MG','Activin + cyclohex','Activin + MG','BMP + cyclohex'};
% % posPerCondition = 4;
% % nWells = 8;

%%
barefname = 'syringeC2C12_7hintervals';
treatmentTime = 7;
posPerCondition = 8;
nWells = 1;

nucChannel = 2;
S4Channel = 1;
           
%% read the MIPs from previous step at time 1 

gridSize = meta.montageGridSize;
pixelOverlap = round(1024*meta.montageOverlap/100);

imgsNuc = {};
imgsS4 = {};

for wellnr = 1:nWells
    conditionPositions = posPerCondition*(wellnr-1)+1:posPerCondition*wellnr;
    if isempty(gridSize) 
        gridSize = [posPerCondition/2 2];
    end

    tmax = meta.nTime;
    for ti = 1:tmax

        disp(['processing time ' num2str(ti)]);

        for pi = conditionPositions

            disp(['reading MIP ' num2str(pi)]);
            % gridSize 1 and 2 may be swapped, I have no way of knowing right now
            [i,j] = ind2sub(gridSize, pi - conditionPositions(1) + 1);

            fname = fullfile(dataDir,'MIP',[barefname sprintf('_MIP_p%.4d_w%.4d.tif',pi-1,nucChannel-1)]);
            imgsNuc{j,i} = double(imread(fname,ti));

            fname = fullfile(dataDir,'MIP',[barefname sprintf('_MIP_p%.4d_w%.4d.tif',pi-1,S4Channel-1)]);
            imgsS4{j,i} = double(imread(fname,ti));
        end

        % stitch together
        if ti == 1 && ~isempty(pixelOverlap)
            % get register positions of upper left corner
            upperleft = registerImageGrid(imgsNuc, pixelOverlap);
        elseif ti == 1 && isempty(pixelOverlap)
            upperleft = {};
            for pi = conditionPositions
                [i,j] = ind2sub(gridSize,pi - conditionPositions(1) + 1);
                upperleft{j,i} = [1+(j-1)*(1024 + 50), 1+(i-1)*(1024 + 50)];
            end
        end
        nucStitched = stitchImageGrid(upperleft, imgsNuc);
        S4Stitched = stitchImageGrid(upperleft, imgsS4);

        % make clean preview (not for quantitative analysis
        nucSmall = imfilter(nucStitched,[1 1]/2);
        nucSmall = nucSmall(1:2:end,1:2:end);
        nucSmall = imadjust(mat2gray(nucSmall));
        nucSmall = uint16((2^16-1)*nucSmall);

        S4Small = imfilter(S4Stitched,[1 1]/2);
        S4Small = S4Small(1:2:end,1:2:end);
        S4Small = imadjust(mat2gray(medfilt2(S4Small,[3 3])));
        S4Small = uint16((2^16-1)*S4Small);

        if ti == 1
            previewS4 = zeros([size(nucSmall) tmax],'uint16');
            previewNuc = zeros([size(nucSmall) tmax],'uint16');
        end
        previewNuc(:,:,ti) = nucSmall;
        previewS4(:,:,ti) = S4Small;
    end

    fname = fullfile(dataDir, ['stichedPreviewNuclei_well' num2str(wellnr) '.tif']);
    imwrite(previewNuc(:,:,1), fname);
    for ti = 2:tmax
        imwrite(previewNuc(:,:,ti), fname,'WriteMode','Append');
    end

    fname = fullfile(dataDir, ['stichedPreviewS4_well' num2str(wellnr) '.tif']);
    imwrite(previewS4(:,:,1), fname);
    for ti = 2:tmax
        imwrite(previewS4(:,:,ti), fname,'WriteMode','Append');
    end
end
% figure, imshow(cat(3,0*nucSmall,S4Small,0*S4Small));
% s = strsplit(meta.timeInterval,' ');
% dt = str2double(s{1});
% unit = s{2};
% t = (ti - treatmentTime)*dt;
% text(100,100,['T = ' num2str(t) s{2}],'Color','white','FontSize',18);

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