clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
dataDir = '/Users/idse/data_tmp/0_kinetics/170610_FRAPcyt';

oibfiles = {'A1002hcyt.oib','A1002hcyt2.oib', 'A1002dot45hcyt.oib','A1003.5h.oib','untreated.oib'};
oibfile = oibfiles{4};
[~,barefname,~] = fileparts(oibfile);

load(fullfile(dataDir,[barefname '_results.mat']));

nucChannel = 2;
S4Channel = 1;

%%
plotFRAPcurves(dataDir, oibfile, results)

%%
% fitting nuclear levels

corrected = false;
fitresults = results;
fitresults.bleachType = 'cytoplasmic';
fitresults.fitType = 'nuclear';

idx = fitresults.FrappedCytoIdx;
refidx = [];

tracesNuc = fitresults.tracesNuc(idx,:);
refTracesNuc = fitresults.tracesNuc(refidx,:);
meanref = mean(refTracesNuc',2);

if corrected 
    tracesNucCorr = bsxfun(@minus,tracesNuc',meanref) + max(meanref);
    minI = min(min(fitresults.tracesCyt'));
    tracesNucNormCorr = (tracesNucCorr - minI)./(max(tracesNucCorr)-minI);
else
    tracesNucCorr = tracesNuc';
    tracesNucNormCorr = fitresults.tracesNucNorm(idx,:)';
end

fitresults.tracesNuc = tracesNucCorr';
fitresults.tracesCyt = fitresults.tracesCyt(idx,:);
fitresults.tracesNucNorm = tracesNucNormCorr';
fitresults.tracesCytNorm = fitresults.tracesCytNorm(idx,:);

[parameters, frapframe, gof] = fitFRAP(fitresults);
fitresults.A = parameters.A;
fitresults.k = parameters.k;
fitresults.B = parameters.B;
fitresults.frapframe = frapframe;
fitresults.gof = gof;

% plot the FRAP curves
plotFRAPcurves(dataDir, oibfile, fitresults)

% FRAP fit
visualizeFRAPfit(fitresults)
%xlim([0 tcut*tres]);
saveas(gcf,fullfile(dataDir, [barefname '_FRAPfitNuc.fig' ]));
saveas(gcf,fullfile(dataDir, [barefname '_FRAPfitNuc.png']));
close;

%%
% fitting cytoplasmic levels

fitresults = results;
fitresults.bleachType = 'cytoplasmic';
fitresults.fitType = 'cytoplasmic';

idx = fitresults.FrappedCytoIdx;
refidx = [1];

tracesCyt = fitresults.tracesCyt(idx,:);
refTracesCyt = fitresults.tracesCyt(refidx,:);
meanref = mean(refTracesCyt',2);
tracesCytCorr = bsxfun(@minus,tracesCyt',meanref) + max(meanref);
tracesNucNormCorr = (tracesCytCorr - min(tracesCytCorr))./(max(tracesCytCorr)-min(tracesCytCorr));

fitresults.tracesNuc = tracesNucCorr';
fitresults.tracesCyt = fitresults.tracesCyt(idx,:);%tracesCytCorr';%
fitresults.tracesNucNorm = tracesNucNormCorr';
fitresults.tracesCytNorm = fitresults.tracesCytNorm(idx,:);%tracesNucNormCorr';%

[parameters, frapframe, gof] = fitFRAP(fitresults);
fitresults.A = parameters.A;
fitresults.k = parameters.k;
fitresults.frapframe = frapframe;
fitresults.gof = gof;

% plot the FRAP curves
%plotFRAPcurves(dataDir, oibfile, fitresults)

% FRAP fit
labels = strread(num2str(find(idx')),'%s');
visualizeFRAPfit(fitresults, labels)
%xlim([0 tcut*tres]);
saveas(gcf,fullfile(dataDir, [barefname '_FRAPfitCyt.fig' ]));
saveas(gcf,fullfile(dataDir, [barefname '_FRAPfitCyt.png']));
close;


