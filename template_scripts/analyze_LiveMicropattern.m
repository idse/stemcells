clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

dataDir = '/Users/idse/data_tmp/170720_MPlive20x';

meta = MetadataAndor(dataDir);

filenameFormat = meta.filename;

% manual metadata
%-------------------

% TODO : modify MetadataAndor to contain all info below

barefname = 'liveMP20x';
treatmentTime = 2;
% workaround to not have to deal with changing stitched previews, actually
% it is 4x4
meta.posPerCondition = 4; 
meta.nWells = 16;

% well 1 (BMP) was treated an hour later after seeing no response
% (pipetting error)
% well 4 (BMP) (13-16) looked morphologically better
% GOOD BMP: 1,3; 13-16
% GOOD BMP+SB: 5-8; 9 (rest have "holes" of flat cells)

nucChannel = 2;
S4Channel = 1;

meta.montageGridSize = [2 2];
meta.montageOverlap = 58; % percentage 

tmax = 167;

% visualize positions
%---------------------

% meta.displayPositions;

% TODO: create merged cellData for montage
% movies of distribution over time

%%
upperleft = stitchedPreviews(dataDir, meta); 
save(fullfile(dataDir,'upperleft'),'upperleft');

%% extract nuclear and cytoplasmic levels

% externally I will have all indices starting at 1
% the Andor offset to start at 0 will be internal

opts = struct(...
                    'dataChannels',     [S4Channel nucChannel],... %'fgChannel',        S4Channel,...
                    'cytoplasmicLevels',true,... 
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'MIPidxDir',        fullfile(dataDir,'MIP'),...
                    'tMax',             5,...
                    'nucShrinkage',     2,...
                    'cytoSize',         5,...
                    'bgMargin',         10,...
                    'NCRcutoff',        [3 Inf]);

opts.cleanupOptions = struct('separateFused', true,...
    'clearBorder',true, 'minAreaStd', 0, 'minSolidity',0.95, 'minArea',100);


%% check that the options are set right

pi = 49;
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
imshow(cat(3, A + s*debugInfo.nucmaskraw, A + s*nucmask, A + s*cytmask));

% %%
% nucmaskraw = debugInfo.segall{2}(:,:,1);
% 
% opts.cleanupOptions = struct('separateFused', true,...
%     'clearBorder',true, 'minAreaStd', 0, 'minSolidity',0.95, 'minArea',100);
% 
% nuclearMask = nucmaskraw;
% options = opts.cleanupOptions;
% options.openSize = 5;
% 
% nuclearMask = bwareaopen(nuclearMask, options.minArea/5);
% nuclearMask = imopen(nuclearMask, strel('disk',options.openSize));
% % fill smaller holes that can appear in nuclear segmentation:
% nuclearMask = ~bwareaopen(~nuclearMask,options.minArea/5); 
% nuclearMaskDefused = separateFusedNuclei(nuclearMask,options);
% 
% figure, imshow(cat(3,double(nucmaskraw),nuclearMask,nuclearMaskDefused))

%% run the analysis on all time points

opts.tMax = tmax;

tic
positions(meta.nPositions) = DynamicPositionAndor();

for pi = 57:meta.nPositions

    positions(pi) = DynamicPositionAndor(meta, pi);
    positions(pi).extractData(dataDir, nucChannel, opts);
    positions(pi).makeTimeTraces();
    save(fullfile(dataDir,'positions'), 'positions');
end
toc

%%
load(fullfile(dataDir,'positions'));
load(fullfile(dataDir,'upperleft'),'upperleft');
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

%% superimpose cell coordinates from different frames

ti = 1;
dmax = 5;
coli = 13;

pixelOverlap = round(meta.xSize*meta.montageOverlap/100);
colstarti = (coli - 1)*4 + 1;
colendi = (coli - 1)*4 + 4;
conditionPositions = colstarti:colendi;

for pi = conditionPositions

    [i,j] = ind2sub(meta.montageGridSize, pi - conditionPositions(1) + 1);
    ci = S4Channel;
    fname = fullfile(dataDir,'MIP',[barefname sprintf(['_MIP_p%.4d_w%.4d.tif'],pi-1,ci-1)]);
    imgs{j,i} = imread(fname,ti);
end
%upperleft = registerImageGrid(imgs, pixelOverlap);
stitched = stitchImageGrid(upperleft{coli}{ti}, imgs);

figure,
imshow(imadjust(stitched),[])
XY = {};
hold on
for pi = conditionPositions
    
    [i,j] = ind2sub(meta.montageGridSize, pi - conditionPositions(1) + 1);
    %XY{pi} = positions(pi).cellData(ti).XY + fliplr(upperleft{ti}{j,i});
    % for old matlab
    XYpi = positions(pi).cellData(ti).XY;
    XY{pi} = XYpi + repmat(fliplr(upperleft{coli}{ti}{j,i}),[size(XYpi,1) 1]);
    scatter(XY{pi}(:,1),XY{pi}(:,2),200,'.');
end
hold off
saveas(gcf, fullfile(dataDir,['superimposed_col' num2str(coli) '.png']));

%%

XYall = {};
for ti = 1:meta.nTime-1
    for pi = conditionPositions
        
        [i,j] = ind2sub(meta.montageGridSize, pi - conditionPositions(1) + 1);
        %XY{pi} = positions(pi).cellData(ti).XY + fliplr(upperleft{ti}{j,i});
        XYpi = positions(pi).cellData(ti).XY;
        XY{pi} = XYpi + repmat(fliplr(upperleft{coli}{ti}{j,i}),[size(XYpi,1) 1]);
    end
    XYall{ti} = cat(1,XY{:});
end

%% make single position object for entire colony

dmax = 5;

radiusMicron = 350;

% objected to store the whole colony
P = DynamicColonyAndor(meta, coli, radiusMicron);

cellDataAll = cat(1,positions(conditionPositions).cellData);
cellDataComb = struct();
    
times = 1:meta.nTime - 1;
for ti = times 

    fprintf('.');
    
    % match cells
    matched = distmat(XYall{ti}) < dmax;
    cells = unique(matched,'rows');
    %cells = cells./sum(cells,2);
    cells = bsxfun(@rdivide, cells, sum(cells,2));
    XYcomb = cells*XYall{ti};

    % remove cells that are not on the pattern
    %XYcm = XYcomb - mean(XYcomb);
    XYcm = bsxfun(@minus, XYcomb, mean(XYcomb));
    bad = sqrt(sum(XYcm.^2,2)) > 1.1*P.radiusPixel;
    cells(bad,:) = [];
    XYcomb = cells*XYall{ti};

    % combine cellData
    cellDataComb(ti).XY  = XYcomb;
    cellDataComb(ti).area = cells*cat(1, cellDataAll(:,ti).area);
    cellDataComb(ti).cytLevel = cells*cat(1, cellDataAll(:,ti).cytLevel);
    cellDataComb(ti).nucLevel = cells*cat(1, cellDataAll(:,ti).nucLevel);
    cellDataComb(ti).background = mean(cat(1, cellDataAll(:,ti).background));
    
    if mod(ti,80)==0
        fprintf('\n');
    end
end
fprintf('\n');

P.cellData = cellDataComb;
P.dataChannels = positions(conditionPositions(1)).dataChannels;
save(fullfile(dataDir,['colony' num2str(coli)]), 'P');

%%
load(fullfile(dataDir,['colony' num2str(coli)]));
cellDataComb = P.cellData;

%%
% load stitched image
times = 160;%[1 10 100 160];
for ti = 1:meta.nTime%times
    for pi = conditionPositions

        [i,j] = ind2sub(meta.montageGridSize, pi - conditionPositions(1) + 1);
        ci = S4Channel;
        fname = fullfile(dataDir,'MIP',[barefname sprintf(['_MIP_p%.4d_w%.4d.tif'],pi-1,ci-1)]);
        imgs{j,i} = imread(fname,ti);
    end
    [stitched, ~] = stitchImageGrid(upperleft{coli}{ti}, imgs);

    % try to exclude some dead cells
    nucI = cellDataComb(ti).nucLevel(:,nucChannel) - cellDataComb(ti).background(nucChannel);
    bad = nucI > 3000; %hist(N,50)

    % calc N:C ratio
    N = cellDataComb(ti).nucLevel(~bad,S4Channel) - cellDataComb(ti).background(S4Channel);
    C = cellDataComb(ti).cytLevel(~bad,S4Channel) - cellDataComb(ti).background(S4Channel);
    R = N./C;
    ncrmin = 0.5;
    ncrmax = 1;
    R(R > ncrmax) = ncrmax;
    R(R < ncrmin) = ncrmin;
    Rnorm = (R - min(R))./(max(R)-min(R));
    colidx = round(Rnorm*255 + 1);

    tic
    %figure,
    imshow(imadjust(stitched),[])
    %saveas(gcf, fullfile(dataDir,['col' num2str(coli) '_t' num2str(ti) '_NCRnooverlay.png']));
    hold on
    colors = jet(256);
    XYcomb = cellDataComb(ti).XY;
    scatter(XYcomb(~bad,1),XYcomb(~bad,2), 800, colors(colidx,:),'.','LineWidth',1.5)
    hold off
    toc 
    %saveas(gcf, fullfile(dataDir,['col' num2str(coli) '_t' num2str(ti) '_NCRoverlay.png']));
    %saveas(gcf, fullfile(dataDir,['col' num2str(coli) '_t' num2str(ti) '_NCRoverlay.pdf']));
    %close;
    axis off
    
    frame{ti} = export_fig(gcf,'-native -m2');
end

%%

[X, Y] = ndgrid(1:size(frame{1},1),1:size(frame{1},2));
R = floor(size(frame{1},1)/2);
mask = (X - R).^2 + (Y-R).^2 > R.^2;
mask3 = repmat(mask,[1 1 3]);
%imshow((X - R).^2 + (Y-R).^2 > R.^2,[])

%%
saveResult = true;
if saveResult
    disp('saving');
    v = VideoWriter(fullfile(dataDir,['signalingPatternX_col' num2str(coli) '.mp4']),'MPEG-4');
    v.FrameRate = 5;
    open(v)
    for ti = 12:167%meta.nTime

        im = frame{ti}(75:656, 85:666,:);
        if ti == 12
            R = floor(size(im,1)/2);
            [X, Y] = ndgrid(1:size(im,1),1:size(im,2));
            mask = (X - R).^2 + (Y-R).^2 > R.^2;
            mask3 = repmat(mask,[1 1 3]);
        end
        im(mask3) = 255;
        %imshow(im)
        writeVideo(v,im)
    end
    close(v);
end

%% make radial profiles

P.makeRadialAvgSeg();

%% visualize radial profiles

figure, 
hold on

times = [20 24 28 32 36 40]*4;
colors = jet(numel(times));
i = 1;

for ti = times

    cytAvg = P.radialProfile(ti).CytAvgSeg(:,S4Channel);
    nucAvg = P.radialProfile(ti).NucAvgSeg(:,S4Channel);
    cytStd = P.radialProfile(ti).CytStdSeg(:,S4Channel);
    nucStd = P.radialProfile(ti).NucStdSeg(:,S4Channel);
    R = P.radialProfile(ti).NucCytRatio(:,S4Channel);
    bins = P.radialProfile(ti).BinEdges(1:end-1)';

%     errorbar(bins, cytAvg, cytStd,'-')
%     hold on
%     errorbar(bins, nucAvg, nucStd,'-r'); 
%     hold off

    plot(R,'Color',colors(i,:),'LineWidth',2)
    i = i + 1;
    %plot(nucAvg./cytAvg);    
end
hold off

%% a little smoothing in time

N = 2;
M = 3;
ratioInTime = P.timeTraces.radAvgRatio{1}(2:167,:)';
ratioInTime = conv2(ratioInTime, ones([M N])/(M*N),'valid');

thr = (1:size(ratioInTime,2))/4;

i = 1;
clf
hold on
for ti = times
    plot(bins(2:end-1)*meta.xres,ratioInTime(:,ti),'Color',colors(i,:),'LineWidth',2);
    i = i+1;
end
hold off
legend(strcat(strread([num2str(round(times/4))],'%s'),' h'),'location','northwest')
axis([0 bins(end-1)*meta.xres 0.6 1])

fs = 24;
xlabel('radial position (micron)', 'FontSize',fs, 'FontWeight','Bold');
ylabel('nuclear:cyto Smad4', 'FontSize',fs, 'FontWeight','Bold');
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')

saveas(gcf, fullfile(dataDir,['col' num2str(coli) '_radialProfiles.png']));

%% plot fixed radial position in time

plot(thr, ratioInTime(1:5:end,:),'LineWidth',2)
xlim([thr(1) thr(end)]);
fs = 20;
xlabel('time', 'FontSize',fs, 'FontWeight','Bold');
ylabel('nuc:cyt', 'FontSize',fs, 'FontWeight','Bold');
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')

legend(strread(num2str(round(bins(2:5:end))'),'%s'),'location','southwest')

saveas(gcf, fullfile(dataDir,['col' num2str(coli) '_differentRadii.png']));

% see derivative of ~ -0.5/60 = -0.008

%% kymograph

figure,
colormap jet
RIT = ratioInTime';
RIT(RIT < ncrmin) = ncrmin;
RIT(RIT > ncrmax) = ncrmax;
RIT = (RIT - 0.5)*2;
imagesc(RIT,'XData',bins*meta.xres,'YData',thr)
axis square
fs = 30;
xlabel('radial position (micron)', 'FontSize',fs, 'FontWeight','Bold');
ylabel('time (hrs)', 'FontSize',fs, 'FontWeight','Bold');
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
h = colorbar;
ylabel(h, 'signaling strength','LineWidth', 2, 'FontSize', fs,'FontWeight', 'bold')
%title('signaling strength');
saveas(gcf, fullfile(dataDir,['col' num2str(coli) '_kymograph.png']));

% 100 microns in 10 hours: 10 micron per hour, one cell per hour

%% time derivatives

dRdt = conv2(ratioInTime', [1 0 -1]'/2,'valid');
imagesc(dRdt','XData',bins*meta.xres,'YData',(1:165)/4)
axis square
plot(thr(2:end),dRdt(:,10))

%% 

figure,
plot(thr, conv2(ratioInTime(20,:),GaussD(2,1)','same'))
axis([thr(1) thr(end) -0.02 0.02]);

%% distributions



%%

% % refine shift between frames based on matching cells
% U = {[0 0]};
% for subpi = 1:3
%     
%     nsubpi = subpi+1;
% 
%     XY1 = XY{subpi};
%     XY2 = XY{nsubpi};
% 
%     [dist, yidx] = pointMatch(XY1, XY2);
%     matchIdx = dist < 10;
%     dXY = XY2(yidx(matchIdx),:) - XY1(matchIdx,:);
%     U{nsubpi} = mean(dXY);
%     XY{nsubpi} = XY{nsubpi} - U{nsubpi};
% end
% 
% % visualize
% clf
% hold on
% for subpi = 1:4
%     scatter(XY{subpi}(:,1), XY{subpi}(:,2))
% end
% axis equal