clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
dataDir = '/Users/idse/data_tmp/0_kinetics/170610_FRAPcyt';

oibfiles = {'A1002hcyt.oib','A1002hcyt2.oib', 'A1002dot45hcyt.oib','A1003.5h.oib','untreated.oib'};

tmaxall = {[],[],[],[],[]};

allresults = {};

descriptions = {'A100 2h', 'A100 2h', 'A100 2.45h','A100 3.5h','untreated'};

nucChannel = 2;
S4Channel = 1;

%% check segmentation parameters

fi = 4;
meta = Metadata(fullfile(dataDir, oibfiles{fi}));

opts = struct(...
                    'dataChannels',     [S4Channel nucChannel],...%'fgChannel',        S4Channel,...
                    'cytoplasmicLevels',true,... 
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'tMax',             meta.nTime-1,...
                    'nucShrinkage',     5,...
                    'cytoSize',         10,...
                    'bgMargin',         10,...
                    'NCRcutoff',        [Inf Inf]);

opts.cleanupOptions = struct('separateFused', true,...
'clearBorder',true, 'minAreaStd', 0, 'minSolidity',0.95, 'minArea',3000);

P = Position(meta.nChannels, oibfiles{fi}, meta.nTime);
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
imshow(cat(3, A, A + s*nucmask, A + s*cytmask));

opts.tMax = meta.nTime;

%% test mask

bla = P.loadSegmentation(fullfile(dataDir,'MIP'), nucChannel);
opts.cleanupOptions = struct('separateFused', true,...
'clearBorder',true, 'minAreaStd', 0, 'minSolidity',0.95, 'minArea',3000);

ti = 50;
mask = bla(:,:,ti)';
imshow(mask)
maskclean = nuclearCleanup(mask, opts.cleanupOptions);
imshow(cat(3,double(mask),maskclean,maskclean))

%% analysis 
% for cytoplasm requires first making mips and segmenting nuclei
% TODO : MAKE PREVIEW OF SEGMENTATION

for fi = 4%:numel(oibfiles)
    
    oibfile = oibfiles{fi};
    [~,barefname,~] = fileparts(oibfile);

    % read the data
    [img, omeMeta] = readStack2(fullfile(dataDir,oibfile));
    
    % extract data
    fname = fullfile(dataDir, oibfile);
    meta = Metadata(fname);
    position = Position(meta.nChannels, fname, meta.nTime-1);
    position.setID(fi);
    
    opts.tMax = meta.nTime-1;
    opts.saveCleanNucMask = true;
    
    position.extractData(dataDir, nucChannel, opts);
    position.makeTimeTraces();
    position.makeTracks();

    % detect FRAP frame
    imgmean = mean(mean(img,1),2);
    S4mean = squeeze(imgmean(:,:,S4Channel,:,1:20));
    [~,frapframe] = min(conv(S4mean, [1 0 -1]','valid'));
    frapframe = frapframe + 2;

    % visualize tracks
    figure, 
    ci = 1; zi = 1;
    imshow(imadjust(mat2gray(img(:,:,ci,zi,frapframe))))
    nTracks = numel(position.timeTraces.trackXY);
    colors = hsv(nTracks);
    hold on
    for i = 1:nTracks
        XY = position.timeTraces.trackXY{i};
        plot(XY(:,1),XY(:,2),'LineWidth',2,'Color',colors(i,:));
        text(XY(1,1) + 30, XY(1,2), num2str(i),'Color', colors(i, :),'FontSize',16);
    end
    axis off
    hold off
    saveas(gcf,fullfile(dataDir,[barefname '_tracks.png']));

    h2 = figure; hold on;
    h1 = figure; hold on;

    colors = hsv(nTracks);
    
    t = meta.timeInterval*(1:meta.nTime);
    FrappedCytoIdx = false([nTracks 1]);
    traces = position.timeTraces.cytLevel;
    legendstr = {};
    
    S4tracesCombinedCyt = [];
    S4tracesCombinedCytNorm = [];
    
    S4tracesCombinedNuc = [];
    S4tracesCombinedNucNorm = [];
    
    for i = 1:nTracks
        
        % identify which cells got the cytoplasm frapped
        trackt = position.timeTraces.trackT{i};
        trackIcyt = position.timeTraces.cytLevel{i}(:,S4Channel);
        trackInuc = position.timeTraces.nucLevel{i}(:,S4Channel);
        
        S4cytTrace = interp1(trackt, trackIcyt, 1:meta.nTime);
        S4nucTrace = interp1(trackt, trackInuc, 1:meta.nTime);
        
        % derivative at frapframe (shift from def of derivative)
        % to determine if the cytoplasm was frapped
        dtrace = conv(S4cytTrace, [1 0 -1]','valid');
        FrappedCytoIdx(i) = dtrace(frapframe-2) < -200;
        legendstr{i} = num2str(i);

        %-----------
        minboth = min(min(S4cytTrace), min(S4nucTrace));

        % cytoplasmic
        S4tracesCombinedCyt = cat(1, S4tracesCombinedCyt, S4cytTrace);
        tracenorm = (S4cytTrace - minboth)/(max(S4cytTrace) - minboth);
        S4tracesCombinedCytNorm = cat(1, S4tracesCombinedCytNorm, tracenorm);

        %nuclear
        S4tracesCombinedNuc = cat(1, S4tracesCombinedNuc, S4nucTrace);
        tracenorm = (S4nucTrace - minboth)/(max(S4nucTrace) - minboth);
        S4tracesCombinedNucNorm = cat(1, S4tracesCombinedNucNorm, tracenorm);
        %-----------

        figure(h1);
        plot(t, S4cytTrace, 'LineWidth',1, 'Color', colors(i,:))
        
        figure(h2);
        plot(t, S4nucTrace, 'LineWidth',1, 'Color', colors(i,:))
    end
    
    figure(h1)
    hold off
    ylabel('cytoplasmic intensity');
    xlabel('time (sec)');
    legend(legendstr);
    saveas(gcf,fullfile(dataDir,[barefname '_tracesCyt.png']));
    
    figure(h2)
    hold off
    ylabel('nuclear intensity');
    xlabel('time (sec)');
    legend(legendstr);
    saveas(gcf,fullfile(dataDir,[barefname '_tracesNuc.png']));
    
    % set up results struct
    results = struct('description', descriptions{fi}, 'tmax', tmaxall(fi));
    results.tres = meta.timeInterval;
    results.xyres = meta.xres;
    % store which are cytofrapped in results
    results.FrappedCytoIdx = FrappedCytoIdx;
    
    % read out profile in the mask
    data = squeeze(img(:,:,S4Channel,zi,:));
    results.tracesCyt = S4tracesCombinedCyt;
    results.tracesCytNorm = S4tracesCombinedCytNorm;
    results.tracesNuc = S4tracesCombinedNuc;
    results.tracesNucNorm = S4tracesCombinedNucNorm;
    
    results.tmax = sum(~isnan(results.tracesCyt),2);
    
    % store results of this video
    allresults{fi} = results;
    save(fullfile(dataDir,[barefname '_results']),'results');
    
    % visualize FRAP regions for diagnostic purposes

    % read FRAP regions
    [x,y,shapeTypes] = readFRAPregions(omeMeta);
    
    % initial state
    colors = lines(numel(x));
    ti = frapframe;
    zi = 1;
    ci = S4Channel;
    imshow(imadjust(mat2gray(img(:,:,ci,zi,ti))))
    hold on
    for shapeIdx = 1:numel(x)
        plot(x{shapeIdx}, y{shapeIdx},'Color',colors(shapeIdx,:),'LineWidth',2)
    end
    hold off
    axis off
    saveas(gcf, fullfile(dataDir, [barefname '_FRAPframe' num2str(ti) '.png']));

    % final state
    zi = 1; 
    IS4 = imadjust(mat2gray(img(:,:,ci,zi,end)));
    if meta.nChannels >1
        Inuc = imadjust(mat2gray(img(:,:,nucChannel,zi,end)));
        InucInit = imadjust(mat2gray(img(:,:,nucChannel,zi,1)));
    else
        Inuc = 0*IS4;
    end
    figure, imshow(cat(3, IS4 + Inuc, IS4 + InucInit, IS4));
    axis off
    saveas(gcf, fullfile(dataDir, [barefname '_FRAPframeFinal.png']));
    close;
end

save(fullfile(dataDir,'results'),'allresults');

