clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
mainDataDir = '/Users/idse/data_tmp/0_kinetics/';

FRAPmetadata

%% copied from FRAP_overview

%short time res: 6 2; 11 4; 
untrIdx = [2 4; 4 1; 6 2; 11 4; 10 5;  18 4; 19 4; 21 10];  % 10 4 is too short
peakIdx = [1 1; 3 1; 3 2; 6 1; 7 1; 12 1; 19 2; 21 1; 21 2]; % 8 2 but some RI mistake
adaptIdx = [3 3; 3 4; 9 2; 7 2; 21 4]; % 19 3 

untrIdxLMB = [2 3; 9 6; 9 7; 11 2; 11 3; 16 3];
peakIdxLMB = [9 1; 10 2; 11 1; 13 1; 14 1; 19 1; 21 5; 21 6];% 2 1, but that one recoveres very little and too fast
adaptIdxLMB = [9 3; 10 1; 16 2; 17 1; 21 8];

allIdx = {untrIdx, peakIdx, adaptIdx, untrIdxLMB, peakIdxLMB, adaptIdxLMB};

%%

% NEXT all fits for complete time interval, see if that is better

for j = 1%:6

for i = 1:size(allIdx{j},1)
% 21 6; 23 4; 1 1; 2 2; 2 3
idx = allIdx{j}(i,:);

% k = 24;
% for i = 6:numel(oibfiles{k})
%     idx = [k i];
    
% from metadata
fi = idx(2);
frapframe = frapframes(idx(1));
dataDir = fullfile(mainDataDir, FRAPdirs{idx(1)});
oibfile = oibfiles{idx(1)}{idx(2)};
tm = tmaxall{idx(1)}{idx(2)};
disp(oibfile);

% reading in FRAP curves

[~,barefname,~] = fileparts(oibfile);
barefname = strrep(barefname,'.','dot');

r = bfGetReader(fullfile(dataDir,oibfile));
omeMeta = r.getMetadataStore();

tcut = r.getSizeT()-1; % cutoff time
tres = double(omeMeta.getPixelsTimeIncrement(0).value);
if j == 5 && i == 4
    % override for a file where ROI are corrupted
    Nfrapped = 2;
else
    Nfrapped = omeMeta.getShapeCount(0)/2;
end
if isempty(tm)
    tm = tcut*ones([1 Nfrapped]);
end

% read the data
channels = 1:2;
img = zeros([r.getSizeY() r.getSizeX() numel(channels) r.getSizeZ() tcut], 'uint16');
for ti = 1:tcut
    for cii = 1:numel(channels)
        for zi = 1:r.getSizeZ()
            img(:,:,cii,zi,ti) = bfGetPlane(r, r.getIndex(zi-1,channels(cii)-1,ti-1)+1);
        end
    end
end
r.close();

clear results
if exist(fullfile(dataDir,'results171217.mat'),'file')
    disp('loading previous results');
    load(fullfile(dataDir,'results171217'),'allresults');
    if numel(allresults) >= fi && ~isempty(allresults{fi})
        results = allresults{fi};
    else
        warning('this one is empty');
    end
end

if ~exist('results','var')
    disp('set up results struct');
    results = struct(   'description', oibfile,...
                        'tmax', {tm},...
                        'frapframe', frapframe);
    results.tres = tres;
    results.xyres = double(omeMeta.getPixelsPhysicalSizeX(0).value);
    results.meanI = squeeze(mean(mean(img(:,:,S4Channel,1,:),1),2));
end

% for seeing how this will affect things 
results.shrink = 0.1;

% read out profile in the mask
data = squeeze(img(:,:,S4Channel,zi,:));
if j == 5 && i == 4
    % override for a file where ROI are corrupted
    results = readFRAPprofile2(data, omeMeta, results, LMB{idx(1)}(idx(2)), [], 2); % override : [], 2);
else
    results = readFRAPprofile2(data, omeMeta, results, LMB{idx(1)}(idx(2))); % override : [], 2);
end

% bleach factor
bg = results.bg;
bleachFactor = medfilt2(results.meanI',[1 5]) - bg;
bleachFactor(1:frapframe) = bleachFactor(frapframe);
bleachFactor = bleachFactor./bleachFactor(frapframe);
results.bleachFactor = bleachFactor;

% cytoplasmic bleach factor 'x' from SI
cytbefore = (mean(results.tracesCyt(:,1:frapframe-1),2) - bg);
cytobleachFactor = (results.tracesCyt - bg)./repmat(cytbefore,[1 size(results.tracesCyt,2)]);
results.cytobleachFactor = mean(cytobleachFactor(:,frapframe:frapframe+10),2);

allresults{fi} = results;

% determine background levels

% Iraw = img(:,:,1,1,1);
% I = mat2gray(Iraw);
% 
% LMB = false;
% 
% if LMB
%     B = I > graythresh(I);
%     I(B) = 0;
%     mask = (I > graythresh(I)) +  B;
% else
%     mask = I > graythresh(I);
% end
% mask = imclose(mask,strel('disk',3));
% mask = imopen(imfill(mask,'holes'), strel('disk',3));
% bgmask = imerode(~mask,strel('disk',30));
% bg = mean(Iraw(bgmask));% std(double(Iraw(bgmask)))]
% 
% imshow(cat(3, mat2gray(Iraw), mat2gray(Iraw) + bgmask, mat2gray(Iraw)))

% visualize FRAP regions for diagnostic purposes

% initiale state
colors = lines(Nfrapped);
zi = 1;
for ti = [frapframe, frapframe-1]
    figure, imshow(imadjust(mat2gray(img(:,:,S4Channel,zi,ti))))
    hold on
    for shapeIdx = 1:Nfrapped
        if strcmp(results.shapeTypes{shapeIdx}, 'Polygon')
            nucp = results.nucxstart{shapeIdx};
            cytp = results.cytxstart{shapeIdx};
            plot(nucp([1:end 1],1), nucp([1:end 1],2),'Color',colors(shapeIdx,:),'LineWidth',2)
            plot(cytp([1:end 1],1), cytp([1:end 1],2),'Color',colors(shapeIdx,:),'LineWidth',2)
        end
    end
    hold off
    saveas(gcf, fullfile(dataDir, [barefname '_FRAPframe' num2str(ti) '.png']));
    close;
end

% final state
ti = tcut;
zi = 1;
figure, imshow(imadjust(mat2gray(img(:,:,S4Channel,zi,ti))));
hold on
for shapeIdx = 1:Nfrapped
    if strcmp(results.shapeTypes{shapeIdx}, 'Polygon')
        nucp = results.nucxend{shapeIdx};
        cytp = results.cytxend{shapeIdx};
        plot(nucp([1:end 1],1), nucp([1:end 1],2),'Color',colors(shapeIdx,:),'LineWidth',2)
        plot(cytp([1:end 1],1), cytp([1:end 1],2),'Color',colors(shapeIdx,:),'LineWidth',2)   
    end
end
hold off
saveas(gcf, fullfile(dataDir, [barefname '_FRAPframeFinal.png']));
close;

newfname = 'results171217.mat';
if exist(fullfile(dataDir,newfname),'file')
    load(fullfile(dataDir,newfname),'allresults');
    allresults{idx(2)} = results;
end

save(fullfile(dataDir,newfname),'allresults');
end
end

%% manual correction for Aafter2Hr2.oif 
% (this has no empty space and ends up with wrong background)
j = 2;
i = 4;
idx = allIdx{j}(i,:);

% from metadata
dataDir = fullfile(mainDataDir, FRAPdirs{idx(1)});
oibfile = oibfiles{idx(1)}{idx(2)};
newfname = 'results171217.mat';

S = load(fullfile(dataDir,newfname),'allresults');
results = S.allresults{idx(2)};
bg = 150;
results.bgempty = bg;
results.tracesNucNorm = (results.tracesNuc - bg)/(mean(results.tracesNuc(1:2)) - bg);
results.tracesCytNorm = (results.tracesCyt - bg)/(mean(results.tracesCyt(1:2)) - bg);
S.allresults{idx(2)} = results;
allresults = S.allresults;
save(fullfile(dataDir,newfname),'allresults');

% need to rerun trace extraction

%% manual correction for Untreated bleach8 2nucl.oib 
% (this has no empty space and ends up with wrong background)

for i = 6:8
    idx = [24 i];

    % from metadata
    dataDir = fullfile(mainDataDir, FRAPdirs{idx(1)});
    oibfile = oibfiles{idx(1)}{idx(2)};
    newfname = 'results171217.mat';

    S = load(fullfile(dataDir,newfname),'allresults');
    results = S.allresults{idx(2)};
    bg = 180;
    results.bgempty = bg;
    results.tracesNucNorm = (results.tracesNuc - bg)/(mean(results.tracesNuc(1:2)) - bg);
    results.tracesCytNorm = (results.tracesCyt - bg)/(mean(results.tracesCyt(1:2)) - bg);
    S.allresults{idx(2)} = results;
    allresults = S.allresults;
    save(fullfile(dataDir,newfname),'allresults');
end

%%
for j = 1:6
for i = 1:size(allIdx{j},1)

    idx = allIdx{j}(i,:);
    disp([num2str(j) ' ' num2str(i) ': ' oibfiles{idx(1)}{idx(2)}])
end
end

%%
%------------------------------------------------------------------------
% fitting
%------------------------------------------------------------------------

pf = '171217';

for j = 1:6

for i = 1:size(allIdx{j},1)
    % 21 6; 23 4; 1 1; 2 2; 2 3
    idx = allIdx{j}(i,:);

% k = 24;
% for i = 1:numel(oibfiles{k})
%     idx = [k i];
    
    % from metadata
    fi = idx(2);
    frapframe = frapframes(idx(1));
    dataDir = fullfile(mainDataDir, FRAPdirs{idx(1)});
    oibfile = oibfiles{idx(1)}{idx(2)};
    tm = tmaxall{idx(1)}{idx(2)};
    load(fullfile(dataDir,['results' pf]),'allresults');

    disp('---------------------------------------------------------');
    disp(oibfile);

    results = allresults{fi};
    %results.tres = 10;
    tmalt = zeros([1 size(results.tracesNuc,1)]) + size(results.tracesNuc,2);
    if j == 1
        tmalt = zeros([1 size(results.tracesNuc,1)]) + round(900/results.tres);
    end
    if ~isempty(tm)
        tm = min(tm, tmalt);
    else
        tm = tmalt;
    end
    results.tmax = tm;

    % redefine normalized curves relative to background
    T = results.tracesCyt'-results.bgempty;
    %T(T<0) = 0;
    T=T./max(T);
    results.tracesCytNorm = T';
    T = results.tracesNuc'-results.bgempty;
    %T(T<0) = 0;
    T=T./max(T);
    results.tracesNucNorm = T';

    % fit nuclear recovery
    results.bleachType = 'nuclear';
    results.fitType = 'nuclear';
    [parameters, ~, gof] = fitFRAP3(results);
    results.A = parameters.A;
    results.B = parameters.B;
    results.k = parameters.k;
    results.frapframe = frapframe;

    disp('bleach corrected----------------------------------------');
    resultsb = results;
    resultsb.tracesNucNorm = resultsb.tracesNucNorm./results.bleachFactor;
    [parameters, ~, gof] = fitFRAP3(resultsb);
    results.Ab = parameters.A;
    results.Bb = parameters.B;
    results.kb = parameters.k;
    results.tracesNucBC = resultsb.tracesNuc;
    results.tracesNucNormBC = resultsb.tracesNucNorm;
    results.gofn = gof;

%     disp('fit cytoplasmic levels');
%     resultsb  = results;
%     resultsb.tracesCytNormBC = results.tracesCytNorm./results.bleachFactor;
%     resultsb.tracesCytNorm = resultsb.tracesCytNormBC;
%     resultsb.fitType = 'cytoplasmic';
%     [parameters, ~, gof] = fitFRAP(resultsb);
%     results.Acb = parameters.A;
%     results.kcb = parameters.k;
%     results.Bcb = parameters.B;
%     results.gofc = gof;

    % store results of this video
    fname = ['results' pf 'c' '.mat'];
    if exist(fullfile(dataDir,fname),'file')
        load(fullfile(dataDir,fname),'allresults');
    end
    allresults{idx(2)} = results;
    save(fullfile(dataDir,fname),'allresults');
end
end 

%%
%------------------------------------------------------------------------
% visualize results
%------------------------------------------------------------------------

pf = '171217c';

resultsDir = fullfile(mainDataDir,'results');
subDirs = fullfile(resultsDir,{'untreated','peak','adapted',...
            'untreated_LMB','peak_LMB','adapted_LMB'});
for j = 1:numel(subDirs)
    if ~exist(subDirs{j},'dir')
        mkdir(subDirs{j});
    end
end

for j = 1:6

for i = 1:size(allIdx{j},1)
    idx = allIdx{j}(i,:);

% k = 24;
% for i = 8:numel(oibfiles{k})
%     idx = [k i]; 
%     subDirs{1} = fullfile(mainDataDir, FRAPdirs{idx(1)});
    
    % from metadata
    fi = idx(2);
    frapframe = frapframes(idx(1));
    dataDir = fullfile(mainDataDir, FRAPdirs{idx(1)});
    oibfile = oibfiles{idx(1)}{idx(2)};
    load(fullfile(dataDir,['results' pf]),'allresults');
    results = allresults{fi};
    tm = results.tmax;
    disp(oibfile);

    traces = allresults{fi}.tracesNuc;
    tracesnorm = allresults{fi}.tracesNucNorm;
    tracesnormBC = allresults{fi}.tracesNucNormBC;
    tcut = size(tracesnorm,2);
    tres = allresults{fi}.tres;

    if isempty(allresults{fi}.tmax)
        Nfrapped = size(allresults{fi}.tracesNuc,1);
        allresults{fi}.tmax = tcut*ones([1 Nfrapped]);
    end

    % filename
    [~,barefname,~] = fileparts(oibfile);
    barefname = strrep(barefname,'.','dot');
    Nfrapped = size(allresults{fi}.tracesNuc,1);
    t = repmat((1:tcut)*tres,[Nfrapped 1]);
    pref = pf;%'171214';
        
%     % FRAP curves
%     figure,
%     plot(t' ,traces');
%     hold on 
%     plot(t' ,0*allresults{fi}.tracesNuc'+allresults{fi}.bgempty);
%     hold off
%     xlabel('time (sec)');
%     ylabel('intensity')
%     %saveas(gcf,fullfile(dataDir, ['FRAPcurvesRaw_' barefname]));
%     saveas(gcf,fullfile(subDirs{j}, [pref '_FRAPcurvesRaw_' FRAPdirs{idx(1)}(1:6) '_' barefname '.png']));
%     
%     plot(t' ,allresults{fi}.tracesCyt');
%     hold on 
%     plot(t' ,0*allresults{fi}.tracesCyt'+allresults{fi}.bgempty);
%     hold off
%     xlabel('time (sec)');
%     ylabel('intensity')
%     %saveas(gcf,fullfile(dataDir, ['FRAPcurvesRaw_' barefname]));
%     saveas(gcf,fullfile(subDirs{j}, [pref '_FRAPcurvesCytRaw_' FRAPdirs{idx(1)}(1:6) '_' barefname '.png']));
%     close;
%     
%     plot(t' ,tracesnorm');
%     xlabel('time (sec)');
%     ylabel('normalized intensity')
%     %saveas(gcf,fullfile(dataDir, ['FRAPcurvesNorm_' barefname]));
%     saveas(gcf,fullfile(subDirs{j}, [pref '_FRAPcurvesNorm_' FRAPdirs{idx(1)}(1:6) '_' barefname '.png']));
%     close;
%     
%     plot(t' ,allresults{fi}.tracesCytNorm');
%     xlabel('time (sec)');
%     ylabel('normalized intensity')
%     %saveas(gcf,fullfile(dataDir, ['FRAPcurvesNorm_' barefname]));
%     saveas(gcf,fullfile(subDirs{j}, [pref '_FRAPcurvesCytNorm_' FRAPdirs{idx(1)}(1:6) '_' barefname '.png']));
%     close;
%     
%     plot(t' ,tracesnormBC');
%     xlabel('time (sec)');
%     ylabel('normalized intensity')
%     %saveas(gcf,fullfile(dataDir, ['FRAPcurvesNorm_' barefname]));
%     saveas(gcf,fullfile(subDirs{j}, [pref '_FRAPcurvesNormBC_' FRAPdirs{idx(1)}(1:6) '_' barefname '.png']));
%     close;
%     %     
%     %     % bleach curve
%     %     figure,
%     %     plot(t' ,allresults{fi}.meanI');
%     %     xlabel('time (sec)');
%     %     ylabel('intensity')
%     %     %saveas(gcf,fullfile(dataDir, ['FRAPcurvesRaw_' barefname]));
%     %     saveas(gcf,fullfile(dataDir, ['BleachCurve_' barefname '.png']));
%     %     close;
% 
%     % bleach factor
%     figure,
%     plot(t' ,allresults{fi}.bleachFactor');
%     xlabel('time (sec)');
%     ylabel('intensity')
%     %saveas(gcf,fullfile(dataDir, ['FRAPcurvesRaw_' barefname]));
%     saveas(gcf,fullfile(subDirs{j}, ['BleachFactor_' FRAPdirs{idx(1)}(1:6) '_' barefname '.png']));
%     close;

    % FRAP fit
    clf
    visualizeFRAPfit2(allresults{fi},[],true)
    %xlim([0 min(size(results.meanI,1)*tres, max(results.tres*tm)*2)])
    %saveas(gcf,fullfile(dataDir, ['FRAPfit_' barefname]));
    saveas(gcf,fullfile(subDirs{j}, [pref '_FRAPfitNuc_' FRAPdirs{idx(1)}(1:6) '_' barefname '.png']));
    close;

    resultsBC = allresults{fi};
    resultsBC.A = resultsBC.Ab;
    resultsBC.B = resultsBC.Bb;
    resultsBC.k = resultsBC.kb;
    resultsBC.tracesNucNorm = resultsBC.tracesNucNormBC;

    visualizeFRAPfit2(resultsBC,[],true)
    %saveas(gcf,fullfile(dataDir, ['FRAPfitNucBC_' barefname]));
    saveas(gcf,fullfile(subDirs{j}, [pref '_FRAPfitNucBC_' FRAPdirs{idx(1)}(1:6) '_' barefname '.png']));
    close;

%     % cytoplasmic fit
%     resultsCyt = allresults{fi};
%     resultsCyt.A = resultsCyt.Acb;
%     resultsCyt.B = resultsCyt.Bcb;
%     resultsCyt.k = resultsCyt.kcb;
%     resultsCyt.fitType = 'cytoplasmic';
%     clf
%     visualizeFRAPfit(resultsCyt)
%     saveas(gcf,fullfile(dataDir, ['FRAPfitCytBC_' barefname]));
%     saveas(gcf,fullfile(dataDir, ['FRAPfitCytBC_' barefname '.png']));
%     close;
end
end
