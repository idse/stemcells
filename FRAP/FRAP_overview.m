clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

dataDir = '/Users/idse/data_tmp/0_kinetics/';

FRAPdirs = {'170209_FRAP','170302_FRAPagain','170510_FRAPafterA',...
            '170519_FRAPctrl','170524_FRAP','170525_FRAP','170529_FRAPAvB',...
            '170530_FRAP','170817'};%,'170605_Frap'};

%FRAPdirs = FRAPdirs(6);

%% load FRAP data from different experiments

results = {};

for i = 1:numel(FRAPdirs)
	
    S = load(fullfile(dataDir,FRAPdirs{i},'results2.mat'));
    results{i} = S.allresults;
    
%     % hack because I'm loading old files, fix
%     for j = 1:numel(results{i})
%         if ~isempty(results{i}{j})
%             results{i}{j}.tracesNuc = results{i}{j}.traces;
%             results{i}{j}.tracesNucNorm = results{i}{j}.tracesnorm;
%             results{i}{j}.bleachType = 'nuclear';
%         end
%     end
end

save(fullfile(dataDir, 'FRAPresults.mat'), 'results');

% %%
% 
% results{3}{3}.good = [1 1];
% results{3}{4}.good = [1 1];
% results{2}{3}.good = [1 0];
% 
% for i = 1:numel(FRAPdirs)
% 	allresults = results{i};
%     save(fullfile(dataDir,FRAPdirs{i},'results2.mat'),'allresults');
% end
% 

%% combine similar condition measurements

load(fullfile(dataDir, 'FRAPresults.mat'));

untrIdx = [6 2; 2 4];
peakIdx = [1 1; 6 1];
adaptIdx = [3 3; 3 4; 9 2; 9 4; 9 5];

untrIdxLMB = [2 3; 9 6; 9 7];
peakIdxLMB = [9 1];
adaptIdxLMB = [9 3];

instruct = struct(  'A',[],'errA',[],'k',[],'errk',[],...
                    'Nmovies',0, 'Ncells',0,...
                    'rinit',[],'rfin',[]);

untreated_in = instruct;
peak_in = instruct;
adapted_in = instruct;

untreatedLMB_in = instruct;
peakLMB_in = instruct;
adaptedLMB_in = instruct;

disp('reading untreated values');
for i = 1:size(untrIdx,1)
    
    result = results{untrIdx(i,1)}{untrIdx(i,2)};
    goodidx = logical(result.good);
    xr = getNCR(result);
    
    untreated_in.Nmovies = untreated_in.Nmovies + 1;
    untreated_in.Ncells = untreated_in.Ncells + sum(goodidx);

    untreated_in.rinit = cat(1, untreated_in.rinit, xr(:,2));
    untreated_in.rfin = cat(1, untreated_in.rfin, xr(:,3));
    untreated_in.A = cat(1, untreated_in.A, result.A(goodidx,1));
    untreated_in.errA = cat(1, untreated_in.errA, result.A(goodidx,1)-result.A(goodidx,2));
    untreated_in.k = cat(1, untreated_in.k, result.k(goodidx,1));
    untreated_in.errk = cat(1, untreated_in.errk, result.k(goodidx,1)-result.k(goodidx,2));
end

disp('reading peak values');
for i = 1:size(peakIdx,1)
    
    result = results{peakIdx(i,1)}{peakIdx(i,2)};
    goodidx = logical(result.good);
    xr = getNCR(result);
    
    peak_in.Nmovies = peak_in.Nmovies + 1;
    peak_in.Ncells = peak_in.Ncells + sum(goodidx);
    
    peak_in.rinit = cat(1, peak_in.rinit, xr(:,2));
    peak_in.rfin = cat(1, peak_in.rfin, xr(:,3));
	peak_in.A = cat(1, peak_in.A, result.A(goodidx,1));
    peak_in.errA = cat(1, peak_in.errA, result.A(goodidx,1)-result.A(goodidx,2));
    peak_in.k = cat(1, peak_in.k, result.k(goodidx,1));
    peak_in.errk = cat(1, peak_in.errk, result.k(goodidx,1)-result.k(goodidx,2));
end

disp('reading adapted values');
for i = 1:size(adaptIdx,1)
    
    result = results{adaptIdx(i,1)}{adaptIdx(i,2)};
    goodidx = logical(result.good);
    xr = getNCR(result);
    
    adapted_in.Nmovies = adapted_in.Nmovies + 1;
    adapted_in.Ncells = adapted_in.Ncells + sum(goodidx);
    
    adapted_in.rinit = cat(1, adapted_in.rinit, xr(:,2));
    adapted_in.rfin = cat(1, adapted_in.rfin, xr(:,3));
	adapted_in.A = cat(1, adapted_in.A, result.A(goodidx,1));
    adapted_in.errA = cat(1, adapted_in.errA, result.A(goodidx,1)-result.A(goodidx,2));
    adapted_in.k = cat(1, adapted_in.k, result.k(goodidx,1));
    adapted_in.errk = cat(1, adapted_in.errk, result.k(goodidx,1)-result.k(goodidx,2));
end
% ---------------LMB------------

disp('reading untreated LMB values');
for i = 1:size(untrIdxLMB,1)
    
    result = results{untrIdxLMB(i,1)}{untrIdxLMB(i,2)};
    goodidx = logical(result.good);
    xr = getNCR(result);
    
    untreatedLMB_in.Nmovies = untreatedLMB_in.Nmovies + 1;
    untreatedLMB_in.Ncells = untreatedLMB_in.Ncells + sum(goodidx);
    
    untreatedLMB_in.rinit = cat(1, untreatedLMB_in.rinit, xr(:,2));
    untreatedLMB_in.rfin = cat(1, untreatedLMB_in.rfin, xr(:,3));
	untreatedLMB_in.A = cat(1, untreatedLMB_in.A, result.A(goodidx,1));
    untreatedLMB_in.errA = cat(1, untreatedLMB_in.errA, result.A(goodidx,1)-result.A(goodidx,2));
    untreatedLMB_in.k = cat(1, untreatedLMB_in.k, result.k(goodidx,1));
    untreatedLMB_in.errk = cat(1, untreatedLMB_in.errk, result.k(goodidx,1)-result.k(goodidx,2));
end

disp('reading LMB peak values');
for i = 1:size(peakIdxLMB,1)
    
    result = results{peakIdxLMB(i,1)}{peakIdxLMB(i,2)};
    goodidx = logical(result.good);
    xr = getNCR(result);
    
    peakLMB_in.Nmovies = peakLMB_in.Nmovies + 1;
    peakLMB_in.Ncells = peakLMB_in.Ncells + sum(goodidx);
    
    peakLMB_in.rinit = cat(1, peakLMB_in.rinit, xr(:,2));
    peakLMB_in.rfin = cat(1, peakLMB_in.rfin, xr(:,3));
	peakLMB_in.A = cat(1, peakLMB_in.A, result.A(goodidx,1));
    peakLMB_in.errA = cat(1, peakLMB_in.errA, result.A(goodidx,1)-result.A(goodidx,2));
    peakLMB_in.k = cat(1, peakLMB_in.k, result.k(goodidx,1));
    peakLMB_in.errk = cat(1, peakLMB_in.errk, result.k(goodidx,1)-result.k(goodidx,2));
end

disp('reading adapted LMB values');
for i = 1:size(adaptIdxLMB,1)
    
    result = results{adaptIdxLMB(i,1)}{adaptIdxLMB(i,2)};
    goodidx = logical(result.good);
    xr = getNCR(result);
    
    adaptedLMB_in.Nmovies = adaptedLMB_in.Nmovies + 1;
    adaptedLMB_in.Ncells = adaptedLMB_in.Ncells + sum(goodidx);
	
    adaptedLMB_in.rinit = cat(1, adaptedLMB_in.rinit, xr(:,2));
    adaptedLMB_in.rfin = cat(1, adaptedLMB_in.rfin, xr(:,3));
    adaptedLMB_in.A = cat(1, adaptedLMB_in.A, result.A(goodidx,1));
    adaptedLMB_in.errA = cat(1, adaptedLMB_in.errA, result.A(goodidx,1)-result.A(goodidx,2));
    adaptedLMB_in.k = cat(1, adaptedLMB_in.k, result.k(goodidx,1));
    adaptedLMB_in.errk = cat(1, adaptedLMB_in.errk, result.k(goodidx,1)-result.k(goodidx,2));
end

disp('display measured parameters');

measured = {'A','k','rinit','rfin'}; % w and w/o LMB
legloc = {'NorthWest','NorthEast','NorthWest','NorthWest'};
yranges = {[],[],[0 20],[0 5]};

for j = 1:2%1:3

    name = measured{j};
    colors = lines(2);
    fs = 30;

    vals = [mean(untreated_in.(name)) mean(peak_in.(name)) mean(adapted_in.(name))]; 
    errorvals = [   std(untreated_in.(name))/sqrt(untreated_in.Ncells-1)...
                    std(peak_in.(name))/sqrt(peak_in.Ncells-1)...
                    std(adapted_in.(name))/sqrt(adapted_in.Ncells-1)]; 
                %[std(untreated_in.(name)) std(peak_in.(name)) std(adapted_in.(name))]; 

    vals2 = [mean(untreatedLMB_in.(name)) mean(peakLMB_in.(name)) mean(adaptedLMB_in.(name))]; 
    errorvals2 = [  std(untreatedLMB_in.(name))/sqrt(untreatedLMB_in.Ncells-1)...
                    std(peakLMB_in.(name))/sqrt(peakLMB_in.Ncells-1)...
                    std(adaptedLMB_in.(name))/sqrt(adaptedLMB_in.Ncells-1)]; 
                %[std(untreatedLMB_in.(name)) std(peakLMB_in.(name)) std(adaptedLMB_in.(name))]; 

    clf
    b = bar(cat(1,vals,vals2)',0.9);%,'FaceColor',[0.1 0.5 0.1]);
    b(1).FaceColor = colors(1,:);
    b(2).FaceColor = colors(2,:);
    hold on
    errorbar([1 2 3]-0.15,vals, errorvals,'k','LineStyle','none','linewidth',3);
    errorbar([1 2 3]+0.15,vals2, errorvals2,'k','LineStyle','none','linewidth',3);
    hold off
    set(gcf,'color','w');
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    box off
    xlim([0.5 3.5]);
    if ~isempty(yranges{j})
        ylim(yranges{j});
    end
    axis square
    xticks([1 2 3]);
    set(gca,'XTickLabel', {'ctrl', 'peak', 'adapt'});
    legend({name, [name '+LMB']},'Location',legloc{j})
    saveas(gcf, fullfile(dataDir,[name '_bargraph.png']));
end

%%
% compare rinit and rfin

colors = lines(2);
fs = 30;
LMB = true;

if ~LMB 
    titlestr = 'no LMB';
    yrange = [0 1.3];
    
    name = 'rinit';
    vals = [mean(untreated_in.(name)) mean(peak_in.(name)) mean(adapted_in.(name))]; 
    errorvals = [   std(untreated_in.(name))/sqrt(untreated_in.Ncells-1)...
                    std(peak_in.(name))/sqrt(peak_in.Ncells-1)...
                    std(adapted_in.(name))/sqrt(adapted_in.Ncells-1)]; 
    name = 'rfin';
    vals2 = [mean(untreated_in.(name)) mean(peak_in.(name)) mean(adapted_in.(name))]; 
    errorvals2 = [  std(untreated_in.(name))/sqrt(untreated_in.Ncells-1)...
                    std(peak_in.(name))/sqrt(peak_in.Ncells-1)... 
                    std(adapted_in.(name))/sqrt(adapted_in.Ncells-1)]; 
else
    titlestr = 'LMB';
    yrange = [0 20];
    
    name = 'rinit';
    vals = [mean(untreatedLMB_in.(name)) mean(peakLMB_in.(name)) mean(adaptedLMB_in.(name))]; 
    errorvals = [   std(untreatedLMB_in.(name))/sqrt(untreatedLMB_in.Ncells-1)...
                    std(peakLMB_in.(name))/sqrt(peakLMB_in.Ncells-1)...
                    std(adaptedLMB_in.(name))/sqrt(adaptedLMB_in.Ncells-1)]; 
    
    name = 'rfin';
    vals2 = [mean(untreatedLMB_in.(name)) mean(peakLMB_in.(name)) mean(adaptedLMB_in.(name))]; 
    errorvals2 = [  std(untreatedLMB_in.(name))/sqrt(untreatedLMB_in.Ncells-1)...
                    std(peakLMB_in.(name))/sqrt(peakLMB_in.Ncells-1)...
                    std(adaptedLMB_in.(name))/sqrt(adaptedLMB_in.Ncells-1)]; 
end

clf
b = bar(cat(1,vals,vals2)',0.9);%,'FaceColor',[0.1 0.5 0.1]);
b(1).FaceColor = colors(1,:);
b(2).FaceColor = colors(2,:);
hold on
bins = (1:3);
errorbar(bins-0.15, vals, errorvals,'k','LineStyle','none','linewidth',3);
errorbar(bins+0.15, vals2, errorvals2,'k','LineStyle','none','linewidth',3);
hold off
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
box off
xlim([0.5 max(bins)+0.5]);
ylim(yrange);
axis square
xticks([1 2 3]);
set(gca,'XTickLabel', {'ctrl', 'peak', 'adapt'});
legend({'rinit', 'rfin'},'Location',legloc{j})
title(titlestr);
saveas(gcf, fullfile(dataDir,['rcompare' titlestr '.png']));

%% display amounts of input data

Nmovies = [untreated_in.Nmovies peak_in.Nmovies adapted_in.Nmovies]
NmoviesLMB = [untreatedLMB_in.Nmovies peakLMB_in.Nmovies adaptedLMB_in.Nmovies]

Ncells = [untreated_in.Ncells peak_in.Ncells adapted_in.Ncells]
NcellsLMB = [untreatedLMB_in.Ncells peakLMB_in.Ncells adaptedLMB_in.Ncells]

%% overlay some FRAP curves

ei = 1;

colors = lines(3);
lw = 4;
shapeIdx = 1;
clf 
hold on

idx = {[1,1], [5,3], [3, 3]};
shift = [1 28 0];

for i = 1:numel(idx)
    
    ei = idx{i}(1);
    fi = idx{i}(2);
    
    if isfield(results{ei}{fi},'tracesnorm')
        tracesnorm = results{ei}{fi}.tracesnorm(shapeIdx,:);
    else
        tracesnorm = results{ei}{fi}.tracesNucNorm(shapeIdx,:);
    end
    tcut = size(tracesnorm,2);
    tres = results{ei}{fi}.tres;
    
    t = ((1:tcut) + shift(i))*tres;
    tp = (1:tcut)*tres;
    
    func = @(p1,p2,x) p1*(1-exp(-x*p2));
    A = results{ei}{fi}.A(shapeIdx,1);
    k = results{ei}{fi}.k(shapeIdx,1);
    tmax = results{ei}{fi}.tmax;
    frapframe = results{ei}{fi}.frapframe;
    plot(t(frapframe:tmax(shapeIdx)),func(A,k,tp(frapframe:tmax(shapeIdx))),...
                        'Color', colors(i,:),'LineWidth',lw)
end
legend({'peak signaling', 'untreated', 'adapted'},'Location','SouthEast')

for i = 1:numel(idx)%%:numel(oibfiles)
    
    ei = idx{i}(1);
    fi = idx{i}(2);
    
    if isfield(results{ei}{fi},'tracesnorm')
        tracesnorm = results{ei}{fi}.tracesnorm(shapeIdx,:);
    else
        tracesnorm = results{ei}{fi}.tracesNucNorm(shapeIdx,:);
    end
    tcut = size(tracesnorm,2);
    tres = results{ei}{fi}.tres;
    
    t = ((1:tcut) + shift(i))*tres;
    tp = (1:tcut)*tres;
    
    % FRAP curves
    plot(t' ,tracesnorm','Color',colors(i,:));
end

for i = 1:numel(idx)
    
    ei = idx{i}(1);
    fi = idx{i}(2);
    
    if isfield(results{ei}{fi},'tracesnorm')
        tracesnorm = results{ei}{fi}.tracesnorm(shapeIdx,:);
    else
        tracesnorm = results{ei}{fi}.tracesNucNorm(shapeIdx,:);
    end
    tcut = size(tracesnorm,2);
    tres = results{ei}{fi}.tres;
    
    t = ((1:tcut) + shift(i))*tres;
    tp = (1:tcut)*tres;
    
    func = @(p1,p2,x) p1*(1-exp(-x*p2));
    A = results{ei}{fi}.A(shapeIdx,1);
    k = results{ei}{fi}.k(shapeIdx,1);
    tmax = results{ei}{fi}.tmax;
    frapframe = results{ei}{fi}.frapframe;
    plot(t(frapframe:tmax(shapeIdx)),func(A,k,tp(frapframe:tmax(shapeIdx))),...
                        'Color', colors(i,:),'LineWidth',lw)
end

hold off
fs = 32;
ylabel('relative recovery', 'FontSize',fs, 'FontWeight','Bold');
xlabel('time (sec)', 'FontSize',fs, 'FontWeight','Bold');
xlim([0 1000]);
ylim([0 0.55]);
set(gcf,'color','w');
set(gca, 'LineWidth', lw);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')

saveas(gcf, fullfile(dataDir, 'FRAPcombined.png'));
saveas(gcf, fullfile(dataDir, 'FRAPcombined.fig'));

%% MODEL FIT actual values

% measured :
% recovery amplitudes A
% nuc:cyt ratios R
% inverse recovery time k
input = {};
output = {};

instruct = {untreated_in, peak_in, adapted_in};
instructp = {untreatedLMB_in, peakLMB_in, adaptedLMB_in};

for i = 1:3

    input{i} = struct();
        
    input{i}.A = mean(instruct{i}.A);
    input{i}.Ap = mean(instructp{i}.A);
    input{i}.R = mean(instruct{i}.rinit);
    input{i}.Rp = mean(instructp{i}.rinit);
    input{i}.k = mean(instruct{i}.k);
    input{i}.kp = mean(instructp{i}.k);

    % measured standard errors
    input{i}.sigA = std(instruct{i}.A)/sqrt(instruct{i}.Ncells-1);
    input{i}.sigAp = std(instructp{i}.A)/sqrt(instructp{i}.Ncells-1);
    input{i}.sigR = std(instruct{i}.rinit)/sqrt(instruct{i}.Ncells-1);
    input{i}.sigRp = std(instructp{i}.rinit)/sqrt(instructp{i}.Ncells-1);
    input{i}.sigk = std(instruct{i}.k)/sqrt(instruct{i}.Ncells-1);
    input{i}.sigkp = std(instructp{i}.k)/sqrt(instructp{i}.Ncells-1);
    if i == 3
        input{i}.sigAp = 0.1;
        input{i}.sigRp = 0.1;
        input{i}.sigkp = 0.01;
    end
    
    output{i} = fitKineticModel2(input{i});
end

%% visualize model fit

outputcomb = cat(2,output{:});
inferred = {'kin','kout','cs','ns'};

for i = 1:4
    
    name = inferred{i};

    errname = ['sig' name];
    fs = 26;
    vals = [outputcomb.(name)];
    errorvals = [outputcomb.(errname)];
    bar(vals,'FaceColor',[0.1 0.5 0.1]);
    hold on
    errorbar(vals, errorvals,...
                            'k','linestyle','none','linewidth',2);
    hold off
    set(gcf,'color','w');
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    box off
    title(name)
    set(gca,'XTickLabel', {'untreated', 'peak', 'adapt'});
    
    saveas(gcf, fullfile(dataDir, [name '_values.png']));
end

%% MODEL FIT estimated values

% measured :
% recovery amplitudes A
% nuc:cyt ratios R
% inverse recovery time k
fakeinput.A = 0.3;
fakeinput.Ap = 0.2;
fakeinput.R = 0.5;
fakeinput.Rp = 3;
fakeinput.k = 0.003;
fakeinput.kp = 0.0008;

% measured standard errors
fakeinput.sigA = 0.1;
fakeinput.sigAp = 0.1;
fakeinput.sigR = 0.1;
fakeinput.sigRp = 0.1;
fakeinput.sigk = 0.0001;
fakeinput.sigkp = 0.00001;

% just to run the framework for now
%untreated_in = input;
%peak_in = input;
%adapted_in = input;

% fitKineticModel2 same as fitKineticModel but with kappa -> kout
untreated_out = fitKineticModel2(fakeinput);
peak_out = untreated_out;
adapted_out = untreated_out;

%% visualize model fit

name = 'ns';

errname = ['sig' name];
fs = 26;
vals = [untreated_out.(name) peak_out.(name) adapted_out.(name)];
errorvals = [untreated_out.(errname) peak_out.(errname) adapted_out.(errname)];
bar(vals,'FaceColor',[0.1 0.5 0.1]);
hold on
errorbar(vals, errorvals,...
                        'k','linestyle','none','linewidth',2);
hold off
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
box off
set(gca,'XTickLabel', {'untreated', 'peak', 'adapt'});


%% visualize model fit differently

name = 'ns';

errname = ['sig' name];
fs = 26;
vals = [0.2 0.2 0];
errorvals = [untreated_out.(errname) peak_out.(errname) adapted_out.(errname)];
clf
%bar(vals,'FaceColor',[0.1 0.5 0.1]);
hold on
errorbar(vals, errorvals,...
                        'b','linewidth',3);
vals = [0 0.05 0.3];                    
errorbar(vals, errorvals,...
                        'r','linewidth',3);
hold off
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
box off
axis([0.9 3 0 0.4]);
xticks([1 2 3]);
set(gca,'XTickLabel', {'untreated', 'peak', 'adapt'});
legend({'nuclear', 'cytoplasmic'},'Location','NorthWest')


%% visualize model fit differently still

name = 'ns';
colors = lines(2);
errname = ['sig' name];
fs = 30;
vals = [0.2 0.2 0];
vals2 = [0 0.05 0.3];  
errorvals = [untreated_out.(errname) peak_out.(errname) adapted_out.(errname)];
clf
b = bar(cat(1,vals,vals2)',0.9);%,'FaceColor',[0.1 0.5 0.1]);
b(1).FaceColor = colors(1,:);
b(2).FaceColor = colors(2,:);
hold on
errorbar([1 2 3]-0.15,vals, errorvals,'k','LineStyle','none','linewidth',3);
errorbar([1 2 3]+0.15,vals2, errorvals,'k','LineStyle','none','linewidth',3);

% errorbar(vals, errorvals,...
%                         'r','linewidth',3);
hold off
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
box off
axis([0.5 3.5 0 0.4]);
xticks([1 2 3]);
set(gca,'XTickLabel', {'untreated', 'peak', 'adapt'});
legend({'import', 'export'},'Location','NorthWest')
