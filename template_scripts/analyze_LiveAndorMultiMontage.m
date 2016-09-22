clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

dataDir = '/Volumes/IdseData/160429_MPRI';

meta = MetadataAndor(dataDir);

filenameFormat = meta.filename;

% manual metadata
%-------------------

% TODO : modify MetadataAndor to contain all info below

barefname = 'MPRI';
treatmentTime = 4;
posPerCondition = 2;
nWells = 3;

nucChannel = 2;
S4Channel = 1;

meta.montageGridSize = [3 3];
meta.montageOverlap = 10;

% visualize positions
%---------------------

% meta.displayPositions;

% TODO: create merged cellData for montage
% movies of distribution over time

           
%% read the MIPs from previous step at time 1 

gridSize = meta.montageGridSize;
pixelOverlap = 60;%[];%round(1024*meta.montageOverlap/100);

imgsNuc = {};
imgsS4 = {};

for wellnr = 3%1:nWells
    
    for coli = 1:posPerCondition
    
	colnr = (wellnr-1)*posPerCondition + coli;
        
    posPerMontage = prod(meta.montageGridSize);
    conditionPositions = posPerMontage*(colnr-1)+1:posPerMontage*colnr;
    if isempty(gridSize) 
        gridSize = [posPerMontage/2 2];
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
        if ~isempty(pixelOverlap) %&& ti == 1
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
        S4Small = imfilter(S4Stitched,[1 1]/2);
        S4Small = S4Small(1:2:end,1:2:end);

        if ti == 1
            
            IlimNuc = stretchlim(nucSmall);
            IlimS4 = stretchlim(S4Small);
            
            previewS4 = zeros([size(nucSmall) tmax],'uint16');
            previewNuc = zeros([size(nucSmall) tmax],'uint16');
        end
        
        % adjust contrast 
        nucSmall = imadjust(nucSmall, IlimNuc);
        S4Small = imadjust(medfilt2(S4Small,[3 3]), IlimS4);

        yidx = 1:min(size(previewNuc,1),size(nucSmall,1));
        xidx = 1:min(size(previewNuc,2),size(nucSmall,2));
        previewNuc(yidx,xidx,ti) = nucSmall(yidx,xidx);
        previewS4(yidx,xidx,ti) = S4Small(yidx,xidx);
    end

    fname = fullfile(dataDir, ['stichedPreviewNuclei_well' num2str(wellnr) '_col' num2str(coli) '.tif']);
    imwrite(previewNuc(:,:,1), fname);
    for ti = 2:tmax
        imwrite(previewNuc(:,:,ti), fname,'WriteMode','Append');
    end

    fname = fullfile(dataDir, ['stichedPreviewS4_well' num2str(wellnr) '_col' num2str(coli) '.tif']);
    imwrite(previewS4(:,:,1), fname);
    for ti = 2:tmax
        imwrite(previewS4(:,:,ti), fname,'WriteMode','Append');
    end

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
                    'nucShrinkage',     2,...
                    'MIPidxDir',        fullfile(dataDir,'MIP'));

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

wellsWanted = [1 2 6 7];
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

    conditionPositions = 4*(wellnr-1)+1:4*wellnr;
    
    ttraceCat = cat(1,positions.timeTraces);
    ttraceCat = ttraceCat(conditionPositions);
    nucMean = mean(cat(2,ttraceCat.nucLevelAvg),2);
    cytMean = mean(cat(2,ttraceCat.cytLevelAvg),2);
    bgMean = mean(cat(2,ttraceCat.background),2);
    ratioMean = (nucMean - bgMean)./(cytMean - bgMean);
    
    bad = any(cat(1,positions(conditionPositions).ncells) < minNCells,1);
    ratioMean(bad) = NaN;
    plot(t, ratioMean,'LineWidth',2,'Color',colors(wellidx,:))

    g0 = [0 1 0];
    if max(ratioMean) > 1
        g0 = g0/max(ratioMean);
    end

    fs = 24;
    xlabel(['time (' unit ')'], 'FontSize',fs,'FontWeight','Bold')
    ylabel('nuclear : cytoplasmic Smad4', 'FontSize',fs,'FontWeight','Bold');

    axis([t(1), t(end)+50, 0.3, 2]);
    set(gcf,'color','w');
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')

    frame{ti} = export_fig(gcf,'-native -m2');
end
hold off
%title('comparison');
set(gca,'FontSize', 16)
legend(conditions(wellsWanted));
if saveResult
    export_fig(['timeTrace_multipleConditions.pdf'],'-native -m2');
end