clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
mainDataDir = '/Users/idse/data_tmp/0_kinetics/';

FRAPmetadata
idx = [14 1];

% from metadata
fi = idx(2);
frapframe = frapframes(idx(1));
dataDir = fullfile(mainDataDir, FRAPdirs{idx(1)});
oibfile = oibfiles{idx(1)}{idx(2)};
tm = tmaxall{idx(1)}{idx(2)};

%% reading in FRAP curves

[~,barefname,~] = fileparts(oibfile);
barefname = strrep(barefname,'.','dot');

r = bfGetReader(fullfile(dataDir,oibfile));
omeMeta = r.getMetadataStore();

tcut = r.getSizeT()-1; % cutoff time
tres = double(omeMeta.getPixelsTimeIncrement(0).value);
Nfrapped = 2;%omeMeta.getShapeCount(0)/2;
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

% set up results struct
results = struct(   'description', oibfile,...
                    'tmax', {tm},...
                    'frapframe', frapframe);
results.tres = tres;
results.xyres = double(omeMeta.getPixelsPhysicalSizeX(0).value);
results.meanI = squeeze(mean(mean(img(:,:,S4Channel,1,:),1),2));
%%
% read out profile in the mask
data = squeeze(img(:,:,S4Channel,zi,:));
results = readFRAPprofile2(data, omeMeta, results); % override : [], 2);

% bleach factor
bg = min(results.tracesNuc(:,frapframe));
bleachFactor = medfilt1(results.meanI',5) - bg;
bleachFactor(1:frapframe) = bleachFactor(frapframe);
bleachFactor = bleachFactor./bleachFactor(frapframe);
results.bleachFactor = bleachFactor;

% cytoplasmic bleach factor 'x' from SI
cytbefore = (mean(results.tracesCyt(:,1:frapframe-1),2) - bg);
cytobleachFactor = (results.tracesCyt - bg)./cytbefore;
results.cytobleachFactor = mean(cytobleachFactor(:,frapframe:frapframe+10),2);

allresults{fi} = results;

% visualize FRAP regions for diagnostic purposes

% initiale state
colors = lines(Nfrapped);
ti = frapframe;
zi = 1;
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

load(fullfile(dataDir,'results'),'allresults');
allresults{idx(2)} = results;
save(fullfile(dataDir,'results'),'allresults');

%% fitting

load(fullfile(dataDir,'results'),'allresults');

tm = tm*0 + round(1000/results.tres);

results = allresults{fi};
results.tmax = tm;

results.bleachType = 'nuclear';
results.fitType = 'nuclear';

[parameters, ~, gof] = fitFRAP(results);
results.A = parameters.A;
results.k = parameters.k;
results.frapframe = frapframe;

% bleach corrected
resultsb = results;
resultsb.tracesNucNorm = resultsb.tracesNucNorm./results.bleachFactor;

[parameters, ~, gof] = fitFRAP(resultsb);
results.Ab = parameters.A;
results.kb = parameters.k;
results.tracesNucBC = resultsb.tracesNuc;
results.tracesNucNormBC = resultsb.tracesNucNorm;
results.gofn = gof;
% 
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
allresults{idx(2)} = results;
save(fullfile(dataDir,'results'),'allresults');

%end


%%% visualize results

%for fi = 5%1:numel(oibfiles)

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
%     
%     % FRAP curves
%     figure,
%     plot(t' ,traces');
%     xlabel('time (sec)');
%     ylabel('intensity')
%     %saveas(gcf,fullfile(dataDir, ['FRAPcurvesRaw_' barefname]));
%     saveas(gcf,fullfile(dataDir, ['FRAPcurvesRaw_' barefname '.png']));
% 
%     plot(t' ,tracesnorm');
%     xlabel('time (sec)');
%     ylabel('normalized intensity')
%     %saveas(gcf,fullfile(dataDir, ['FRAPcurvesNorm_' barefname]));
%     saveas(gcf,fullfile(dataDir, ['FRAPcurvesNorm_' barefname '.png']));
%     close;
% 
%     plot(t' ,tracesnormBC');
%     xlabel('time (sec)');
%     ylabel('normalized intensity')
%     %saveas(gcf,fullfile(dataDir, ['FRAPcurvesNorm_' barefname]));
%     saveas(gcf,fullfile(dataDir, ['FRAPcurvesNormBC_' barefname '.png']));
%     close;
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
    saveas(gcf,fullfile(dataDir, ['BleachFactor_' barefname '.png']));
    close;

% FRAP fit
clf
visualizeFRAPfit(allresults{fi})
xlim([0 min(size(results.meanI,1)*tres, max(results.tres*tm)*2)])
%saveas(gcf,fullfile(dataDir, ['FRAPfit_' barefname]));
saveas(gcf,fullfile(dataDir, ['FRAPfitNuc_' barefname '.png']));
%close;

resultsBC = allresults{fi};
resultsBC.A = resultsBC.Ab;
resultsBC.k = resultsBC.kb;
resultsBC.tracesNucNorm = resultsBC.tracesNucNormBC;

%     visualizeFRAPfit(resultsBC)
%     saveas(gcf,fullfile(dataDir, ['FRAPfitNucBC_' barefname]));
%     saveas(gcf,fullfile(dataDir, ['FRAPfitNucBC_' barefname '.png']));
%     close;
%     
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
