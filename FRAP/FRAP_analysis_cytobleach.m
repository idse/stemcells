clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
mainDataDir = '/Users/idse/data_tmp/0_kinetics/';

FRAPmetadata

%% cyto indices

untrIdx = [12 1; 25 4; 26 4; 26 5; 26 6; 27 3; 27 4; 28 1; 29 1; 30 1];
peakIdx = [25 1; 25 2; 26 1; 27 1];
adaptIdx = [25 3; 26 2; 26 3; 27 2];

untrIdxLMB = [14 1];
peakIdxLMB = [];
adaptIdxLMB = [];

allIdx = {untrIdx, peakIdx, adaptIdx, untrIdxLMB, peakIdxLMB, adaptIdxLMB};

%%

% NEXT all fits for complete time interval, see if that is better

j = 3%1:3

for i = 2%4:size(allIdx{j},1)

idx = allIdx{j}(i,:);

% from metadata
fi = idx(2);
frapframe = frapframes(idx(1));
dataDir = fullfile(mainDataDir, FRAPdirs{idx(1)});
oibfile = cytoibfiles{idx(1)}{idx(2)};
tm = tmaxall{idx(1)}{idx(2)};
disp(oibfile);

% reading in FRAP curves

[~,barefname,~] = fileparts(oibfile);
barefname = strrep(barefname,'.','dot');

r = bfGetReader(fullfile(dataDir,oibfile));
omeMeta = r.getMetadataStore();
tcut = r.getSizeT()-1; % cutoff time
tres = double(omeMeta.getPixelsTimeIncrement(0).value);

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
clear allresults
if exist(fullfile(dataDir,'resultsCyt.mat'),'file')
    disp('loading previous results');
    load(fullfile(dataDir,'resultsCyt'),'allresults');
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
results = readFRAPprofileCytbleach(img, results);

% set cutoff time
Nfrapped = size(results.tracesNuc,1);
if isempty(tm)
    tm = tcut*ones([1 Nfrapped]);
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

% visualize FRAP regions for diagnostic purposes

% initiale state
colors = lines(Nfrapped);
zi = 1;
for ti = [frapframe, frapframe-1]
    figure, imshow(imadjust(mat2gray(img(:,:,S4Channel,zi,ti))))
    hold on
    for shapeIdx = 1:Nfrapped

        nucp = results.nucxstart{shapeIdx};
        cytp = results.cytxstart{shapeIdx};
        plot(nucp([1:end 1],1), nucp([1:end 1],2),'Color',colors(shapeIdx,:),'LineWidth',2)
        plot(cytp([1:end 1],1), cytp([1:end 1],2),'Color',colors(shapeIdx,:),'LineWidth',2)
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
    nucp = results.nucxend{shapeIdx};
    cytp = results.cytxend{shapeIdx};
    plot(nucp([1:end 1],1), nucp([1:end 1],2),'Color',colors(shapeIdx,:),'LineWidth',2)
    plot(cytp([1:end 1],1), cytp([1:end 1],2),'Color',colors(shapeIdx,:),'LineWidth',2)   
end
hold off
saveas(gcf, fullfile(dataDir, [barefname '_FRAPframeFinal.png']));
close;

newfname = 'resultsCyt.mat';
if exist(fullfile(dataDir,newfname),'file')
    load(fullfile(dataDir,newfname),'allresults');
    allresults{idx(2)} = results;
end

save(fullfile(dataDir,newfname),'allresults');
end
%end

%%
for j = 1:3
for i = 1:size(allIdx{j},1)

    idx = allIdx{j}(i,:);
    disp([num2str(j) ' ' num2str(i) ': ' cytoibfiles{idx(1)}{idx(2)}])
end
end

%%
%------------------------------------------------------------------------
% fitting
%------------------------------------------------------------------------

pf = '';

for j = 1:3

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
    oibfile = cytoibfiles{idx(1)}{idx(2)};
    tm = tmaxall{idx(1)}{idx(2)};
    load(fullfile(dataDir,['resultsCyt' pf]),'allresults');

    disp('---------------------------------------------------------');
    disp([dataDir '/' oibfile]);

    results = allresults{fi};
    %results.tres = 10;
    tmalt = zeros([1 size(results.tracesNuc,1)]) + size(results.tracesNuc,2);
    if j == 1 
        tmalt = zeros([1 size(results.tracesNuc,1)]) + min(size(results.tracesNuc,2),round(1500/results.tres));
    end
    if ~isempty(tm)
        tm = min(tm, tmalt);
    else
        tm = tmalt;
    end
    results.tmax = tm;

    % redefine normalized curves relative to background
    T = results.tracesCyt'-results.bg;
    %T(T<0) = 0;
    T=T./max(T);
    results.tracesCytNorm = T';
    T = results.tracesNuc'-results.bg;
    %T(T<0) = 0;
    T=T./max(T);
    results.tracesNucNorm = T';

%     % fit nuclear recovery
%     if j==1 && i==1 % an exception
%         results.frapframe = 2;
%     end
    
    results.bleachType = 'cytoplasmic';
    results.fitType = 'nuclear';
    [parameters, ~, gof] = fitFRAP3(results);
    results.A = parameters.A;
    results.B = parameters.B;
    results.k = parameters.k;
    results.frapframe = frapframe;

    disp('bleach corrected----------------------------------------');
    measuredBleach = results.bleachFactor;
    bleachfac = exp(-(1:numel(results.bleachFactor))/750);%results.bleachFactor;
    disp(['total bleach: ' num2str(bleachfac(end)), '  vs measured: ' num2str(measuredBleach(end))]);
    resultsb = results;
    resultsb.tracesNucNorm = resultsb.tracesNucNorm./bleachfac;
    [parameters, ~, gof] = fitFRAP3(resultsb);
    results.Ab = parameters.A;
    results.Bb = parameters.B;
    results.kb = parameters.k;
    results.tracesNucBC = resultsb.tracesNuc;
    results.tracesNucNormBC = resultsb.tracesNucNorm;
    results.gofn = gof;

    disp('fit cytoplasmic levels----------------------------------');
    resultsb  = results;
    resultsb.fitType = 'cytoplasmic';
    [parameters, ~, gof] = fitFRAP3(resultsb);
    results.Ac = parameters.A;
    results.kc = parameters.k;
    results.Bc = parameters.B;
    results.gofc = gof;

    disp('cytoplasmic bleach corrected');
    resultsb  = results;
    resultsb.tracesCytNormBC = results.tracesCytNorm./bleachfac;
    resultsb.tracesCytNorm = resultsb.tracesCytNormBC;
    resultsb.fitType = 'cytoplasmic';
    [parameters, ~, gof] = fitFRAP3(resultsb);
    results.Acb = parameters.A;
    results.kcb = parameters.k;
    results.Bcb = parameters.B;
    results.tracesCytBC = (results.tracesCyt - results.bg)./results.bleachFactor + results.bg;
    results.tracesCytNormBC = resultsb.tracesCytNormBC;
    results.gofc = gof;

    % store results of this video
    fname = ['resultsCyt' pf '.mat'];
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

pf = '';

resultsDir = fullfile(mainDataDir,'resultsCytobleach');
subDirs = fullfile(resultsDir,{'untreated','peak','adapted'});
for j = 1:numel(subDirs)
    if ~exist(subDirs{j},'dir')
        mkdir(subDirs{j});
	end
end

for j = 1:3

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
    oibfile = cytoibfiles{idx(1)}{idx(2)};
    load(fullfile(dataDir,['resultsCyt' pf]),'allresults');
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
        
    % FRAP curves
    figure,
    plot(t' ,traces');
    hold on 
    plot(t' ,0*allresults{fi}.tracesNuc'+allresults{fi}.bgempty);
    hold off
    xlabel('time (sec)');
    ylabel('intensity')
    %saveas(gcf,fullfile(dataDir, ['FRAPcurvesRaw_' barefname]));
    saveas(gcf,fullfile(subDirs{j}, [pref '_FRAPcurvesRaw_' FRAPdirs{idx(1)}(1:6) '_' barefname '.png']));
    
    plot(t' ,allresults{fi}.tracesCyt');
    hold on 
    plot(t' ,0*allresults{fi}.tracesCyt'+allresults{fi}.bgempty);
    hold off
    xlabel('time (sec)');
    ylabel('intensity')
    %saveas(gcf,fullfile(dataDir, ['FRAPcurvesRaw_' barefname]));
    saveas(gcf,fullfile(subDirs{j}, [pref '_FRAPcurvesCytRaw_' FRAPdirs{idx(1)}(1:6) '_' barefname '.png']));
    close;
    
    plot(t' ,tracesnorm');
    xlabel('time (sec)');
    ylabel('normalized intensity')
    %saveas(gcf,fullfile(dataDir, ['FRAPcurvesNorm_' barefname]));
    saveas(gcf,fullfile(subDirs{j}, [pref '_FRAPcurvesNorm_' FRAPdirs{idx(1)}(1:6) '_' barefname '.png']));
    close;
    
    plot(t' ,allresults{fi}.tracesCytNorm');
    xlabel('time (sec)');
    ylabel('normalized intensity')
    %saveas(gcf,fullfile(dataDir, ['FRAPcurvesNorm_' barefname]));
    saveas(gcf,fullfile(subDirs{j}, [pref '_FRAPcurvesCytNorm_' FRAPdirs{idx(1)}(1:6) '_' barefname '.png']));
    close;
    
    plot(t' ,tracesnormBC');
    xlabel('time (sec)');
    ylabel('normalized intensity')
    %saveas(gcf,fullfile(dataDir, ['FRAPcurvesNorm_' barefname]));
    saveas(gcf,fullfile(subDirs{j}, [pref '_FRAPcurvesNormBC_' FRAPdirs{idx(1)}(1:6) '_' barefname '.png']));
    close;
    %     
    %     % bleach curve
    %     figure,
    %     plot(t' ,allresults{fi}.meanI');
    %     xlabel('time (sec)');
    %     ylabel('intensity')
    %     %saveas(gcf,fullfile(dataDir, ['FRAPcurvesRaw_' barefname]));
    %     saveas(gcf,fullfile(dataDir, ['BleachCurve_' barefname '.png']));
    %     close;

    % bleach factor
    figure,
    plot(t' ,allresults{fi}.bleachFactor');
    xlabel('time (sec)');
    ylabel('intensity')
    %saveas(gcf,fullfile(dataDir, ['FRAPcurvesRaw_' barefname]));
    saveas(gcf,fullfile(subDirs{j}, ['BleachFactor_' FRAPdirs{idx(1)}(1:6) '_' barefname '.png']));
    close;

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

    % cytoplasmic fit
    resultsCyt = allresults{fi};
    resultsCyt.A = resultsCyt.Ac;
    resultsCyt.B = resultsCyt.Bc;
    resultsCyt.k = resultsCyt.kc;
    resultsCyt.fitType = 'cytoplasmic';
    clf
    visualizeFRAPfit2(resultsCyt)
    %saveas(gcf,fullfile(dataDir, ['FRAPfitCytBC_' barefname]));
    saveas(gcf,fullfile(subDirs{j}, [pref '_FRAPfitCyt_' FRAPdirs{idx(1)}(1:6) '_' barefname '.png']));

    % cytoplasmic fit BC
    resultsCyt = allresults{fi};
    resultsCyt.A = resultsCyt.Acb;
    resultsCyt.B = resultsCyt.Bcb;
    resultsCyt.k = resultsCyt.kcb;
    resultsCyt.tracesCytNorm = resultsCyt.tracesCytNormBC;
    resultsCyt.fitType = 'cytoplasmic';
    clf
    visualizeFRAPfit2(resultsCyt)
    %saveas(gcf,fullfile(dataDir, ['FRAPfitCytBC_' barefname]));
    saveas(gcf,fullfile(subDirs{j}, [pref '_FRAPfitCytBC_' FRAPdirs{idx(1)}(1:6) '_' barefname '.png']));
    %close;
end
end
