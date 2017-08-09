clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

dataDir = '/Users/idse/data_tmp/0_kinetics/';

FRAPdirs = {'170209_FRAP','170302_FRAPagain','170510_FRAPafterA',...
            '170519_FRAPctrl','170524_FRAP','170525_FRAP','170529_FRAPAvB',...
            '170530_FRAP'};%,'170605_Frap'};

%% load FRAP data from different experiments

results = {};

for i = 1:numel(FRAPdirs)
	
    S = load(fullfile(dataDir,FRAPdirs{i},'results.mat'));
    results{i} = S.allresults;
    
    % hack because I'm loading old files, fix
    for j = 1:numel(results{i})
        if ~isempty(results{i}{j})
            results{i}{j}.tracesNuc = results{i}{j}.traces;
            results{i}{j}.tracesNucNorm = results{i}{j}.tracesnorm;
            results{i}{j}.bleachType = 'nuclear';
        end
    end
end

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
input.Am = 0.3;
input.Apm = 0.2;
input.Rm = 0.5;
input.Rpm = 3;
input.km = 0.003;
input.kpm = 0.0008;

% measured standard errors
input.sigAm = 0.1;
input.sigApm = 0.1;
input.sigRm = 0.1;
input.sigRpm = 0.1;
input.sigkm = 0.01;
input.sigkpm = 0.001;

% just to run the framework for now
untreated_input = input;
peak_input = input;
adapted_input = input;

untreated_output = fitKineticModel(untreated_input);
peak_output = fitKineticModel(untreated_input);
adapted_output = fitKineticModel(untreated_input);

%% visualize model fit

fs = 26;
bar([0.4 0.1 0.1]);
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
box off
set(gca,'XTickLabel', {'untreated', 'peak', 'adapt'})

%%

bar([1/0.1 1/0.4 1/0.4]);
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
box off
set(gca,'XTickLabel', {'before', 'peak', 'adapt'})
