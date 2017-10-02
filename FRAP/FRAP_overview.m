clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

dataDir = '/Users/idse/data_tmp/0_kinetics/';

FRAPdirs = {'170209_FRAP','170302_FRAPagain','170510_FRAPafterA',...
            '170519_FRAPctrl','170524_FRAP','170525_FRAP','170529_FRAPAvB',...
            '170530_FRAP'};%,'170605_Frap'};

FRAPdirs = FRAPdirs(6);

%% load FRAP data from different experiments

results = {};

for i = 1:numel(FRAPdirs)
	
    S = load(fullfile(dataDir,FRAPdirs{i},'results.mat'));
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

%% combine similar condition measurements

untrIdx = [1 2];
peakIdx = [1 1];

instruct = struct('A',[],'errA',[],'k',[],'errk',[]);

untreated_in = instruct;
peak_in = instruct;
adapted_in = instruct;

untreatedLMB_in = instruct;
peakLMB_in = instruct;
adaptedLMB_in = instruct;

for i = 1:size(untrIdx,1)
    
    result = results{untrIdx(1)}{untrIdx(2)};
    goodidx = logical(result.good);
	untreated_in.A = cat(1, untreated_in.A, result.A(goodidx,1));
    untreated_in.errA = cat(1, untreated_in.errA, result.A(goodidx,1)-result.A(goodidx,2));
    untreated_in.k = cat(1, untreated_in.k, result.k(goodidx,1));
    untreated_in.errk = cat(1, untreated_in.errk, result.k(goodidx,1)-result.k(goodidx,2));
end

for i = 1:size(peakIdx,1)
    
    result = results{peakIdx(1)}{peakIdx(2)};
    goodidx = logical(result.good);
	peak_in.A = cat(1, peak_in.A, result.A(goodidx,1));
    peak_in.errA = cat(1, peak_in.errA, result.A(goodidx,1)-result.A(goodidx,2));
    peak_in.k = cat(1, peak_in.k, result.k(goodidx,1));
    peak_in.errk = cat(1, peak_in.errk, result.k(goodidx,1)-result.k(goodidx,2));
end

%% display measured parameters

measured = {'A','k','r'}; % w and w/o LMB

name = 'A';
colors = lines(2);
fs = 30;

vals = [mean(untreated_in.(name)) mean(peak_in.(name)) mean(adapted_in.(name))]; 
errorvals = [std(untreated_in.(name)) std(peak_in.(name)) std(adapted_in.(name))]; 

vals2 = [mean(untreatedLMB_in.(name)) mean(peakLMB_in.(name)) mean(adaptedLMB_in.(name))]; 
errorvals2 = [std(untreatedLMB_in.(name)) std(peakLMB_in.(name)) std(adaptedLMB_in.(name))]; 

clf
b = bar(cat(1,vals,vals2)',0.9);%,'FaceColor',[0.1 0.5 0.1]);
b(1).FaceColor = colors(1,:);
b(2).FaceColor = colors(2,:);
hold on
errorbar([1 2 3]-0.15,vals, errorvals,'k','LineStyle','none','linewidth',3);
%errorbar([1 2 3]+0.15,vals2, errorvals2,'k','LineStyle','none','linewidth',3);
hold off
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
box off
axis([0.5 3.5 0 0.4]);
xticks([1 2 3]);
set(gca,'XTickLabel', {'untreated', 'peak', 'adapt'});
legend({'A', 'A+LMB'},'Location','NorthWest')

%% get alpha -> need cytoplasm at time of FRAP

[results{1}{1}.tracesNuc(1,4) results{1}{1}.tracesNuc(1,end)]
[results{1}{1}.cytstart{1} results{1}{1}.cytend{1}]

ACN = results{1}{1}.cytstart{1} - results{1}{1}.cytend{1};
ANN = results{1}{1}.tracesNuc(1,end) - results{1}{1}.tracesNuc(1,4);

dN/dC

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
    
    tracesnorm = results{ei}{fi}.tracesnorm(shapeIdx,:);
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
    
    tracesnorm = results{ei}{fi}.tracesnorm(shapeIdx,:);
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
    
    tracesnorm = results{ei}{fi}.tracesnorm(shapeIdx,:);
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
fs = 26;
ylabel('nuclear fluorescence recovery', 'FontSize',fs, 'FontWeight','Bold');
xlabel('time (sec)', 'FontSize',fs, 'FontWeight','Bold');
xlim([0 1000]);
ylim([0 0.61]);
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')

saveas(gcf, fullfile(dataDir, 'FRAPcombined.png'));
saveas(gcf, fullfile(dataDir, 'FRAPcombined.fig'));

%% MODEL FIT

% measured :
% recovery amplitudes A
% nuc:cyt ratios R
% inverse recovery time k
input.A = 0.3;
input.Ap = 0.2;
input.R = 0.5;
input.Rp = 3;
input.k = 0.003;
input.kp = 0.0008;

% measured standard errors
input.sigA = 0.1;
input.sigAp = 0.1;
input.sigR = 0.1;
input.sigRp = 0.1;
input.sigk = 0.0001;
input.sigkp = 0.00001;

% just to run the framework for now
untreated_in = input;
peak_in = input;
adapted_in = input;

% fitKineticModel2 same as fitKineticModel but with kappa -> kout
untreated_out = fitKineticModel2(untreated_in);
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
