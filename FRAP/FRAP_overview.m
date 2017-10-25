clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
dataDir = '/Users/idse/data_tmp/0_kinetics/';

FRAPmetadata

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

results{1}{1}.good = [0 1];
results{1}{2}.good = [1];

results{2}{3}.good = [1 0];

results{3}{3}.good = [1 1];
results{3}{4}.good = [1 1];

results{6}{1}.good = [1 1];
results{6}{2}.good = [1 1];

results{7}{1}.good = [1 1];
results{7}{2}.good = [1 1 1 1 1 1];

results{8}{2}.good = [1 1];

results{10}{4}.good = [1 1 1 1 0]; 
results{10}{5}.good = [0 1 1 1]; 

results{12}{1}.good = [0 1];

results{13}{1}.good = [1]; % A50 LMB 2h movie that didn't fully bleach

save(fullfile(dataDir, 'FRAPresults.mat'), 'results');

%%
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

untrIdx = [2 4; 6 2; 10 4; 10 5; 11 4]; % 5 2 if I want to make square work
peakIdx = [1 1; 6 1; 7 1; 12 1]; % 8 2 but some RI mistake
adaptIdx = [3 3; 3 4; 9 2; 9 4; 9 5; 7 2];

untrIdxLMB = [2 3; 9 6; 9 7; 11 3];
peakIdxLMB = [9 1; 10 2; 11 2; 13 1];
adaptIdxLMB = [9 3; 10 1];

allIdx = {untrIdx, peakIdx, adaptIdx, untrIdxLMB, peakIdxLMB, adaptIdxLMB};
allOut = {};
allLabels = {'untreated', 'peak', 'adapted',...
                    'untreatedLMB', 'peakLMB', 'adaptedLMB'};

debugplots = false;

for i = 1:numel(allIdx)
    legendstr = FRAPdirs(allIdx{i}(:,1));
    allOut{i} = makeInStruct(results, allIdx{i}, debugplots, legendstr, dataDir, allLabels{i});
    close all;
end

untreated_in = allOut{1};
peak_in = allOut{2}; 
adapted_in = allOut{3}; 

untreatedLMB_in = allOut{4};
peakLMB_in = allOut{5};
adaptedLMB_in = allOut{6};

%% display amounts of input data

Nmovies = [untreated_in.Nmovies peak_in.Nmovies adapted_in.Nmovies]
NmoviesLMB = [untreatedLMB_in.Nmovies peakLMB_in.Nmovies adaptedLMB_in.Nmovies]

Ncells = [untreated_in.Ncells peak_in.Ncells adapted_in.Ncells]
NcellsLMB = [untreatedLMB_in.Ncells peakLMB_in.Ncells adaptedLMB_in.Ncells]

%%
disp('display measured parameters');

measured = {'A','Acorr','k','x'};% 'rinit','rfin', % w and w/o LMB
legloc = {'NorthWest','NorthEast','NorthWest','NorthWest'};
yranges = {[],[],[]};%[0 20],[0 5]};

for j = 1:3%1:3

    figure, 
    
    name = measured{j};
    colors = lines(2);
    fs = 30;

    vals = {untreated_in.(name), peak_in.(name), adapted_in.(name)}; 
    means = [mean(untreated_in.(name)) mean(peak_in.(name)) mean(adapted_in.(name))]; 
    errorvals = [   std(untreated_in.(name))/sqrt(untreated_in.Ncells-1)...
                    std(peak_in.(name))/sqrt(peak_in.Ncells-1)...
                    std(adapted_in.(name))/sqrt(adapted_in.Ncells-1)]; 
                %[std(untreated_in.(name)) std(peak_in.(name)) std(adapted_in.(name))]; 
    
    vals2 = {untreatedLMB_in.(name), peakLMB_in.(name), adaptedLMB_in.(name)}; 
    means2 = [mean(untreatedLMB_in.(name)) mean(peakLMB_in.(name)) mean(adaptedLMB_in.(name))]; 
    errorvals2 = [  std(untreatedLMB_in.(name))/sqrt(untreatedLMB_in.Ncells-1)...
                    std(peakLMB_in.(name))/sqrt(peakLMB_in.Ncells-1)...
                    std(adaptedLMB_in.(name))/sqrt(adaptedLMB_in.Ncells-1)]; 
                %[std(untreatedLMB_in.(name)) std(peakLMB_in.(name)) std(adaptedLMB_in.(name))]; 

    clf
    b = bar(cat(1,means,means2)',0.9);%,'FaceColor',[0.1 0.5 0.1]);
    b(1).FaceColor = colors(1,:);
    b(2).FaceColor = colors(2,:);
    hold on
    errorbar([1 2 3]-0.15,means, errorvals,'k','LineStyle','none','linewidth',3);
    errorbar([1 2 3]+0.15,means2, errorvals2,'k','LineStyle','none','linewidth',3);
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
    
    
    % box plots
    figure, 
    valsvec = [];
    idxvec = [];
    coloridx = [];
    colormap lines;
    lw = 3;
    for i = 1:3
        valsvec = cat(1,valsvec, vals{i}, vals2{i});
        coloridx = cat(1,coloridx, vals{i}*0+1, vals2{i}*0 + 2);
        idxvec = cat(1,idxvec, vals{i}*0+2*i-1, vals2{i}*0+2*i);
    end 
    h = boxplot(valsvec, idxvec,'notch','on','Widths',0.3,...
                            'ColorGroup',coloridx,'Colors',lines(2));
    set(h,'linew',lw);
    set(gca, 'LineWidth', lw);
    set(gcf,'color','w');
    if ~isempty(yranges{j})
        ylim(yranges{j});
    end
    xticks([1.5 3.5 5.5]);
    set(gca,'XTickLabel', {'ctrl', 'peak', 'adapt'});
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    box off
    axis square
    saveas(gcf, fullfile(dataDir,[name '_boxplot.png']));
end

%%
% compare rinit and rfin

for LMB = false%[true false]
    figure,
    colors = lines(2);
    fs = 30;
    %  LMB = true;

    if ~LMB 
        titlestr = 'no LMB';
        yrange = [0 1];

        name = 'ncrinit';
        vals = {untreated_in.(name) peak_in.(name) adapted_in.(name)};
        means = [mean(untreated_in.(name)) mean(peak_in.(name)) mean(adapted_in.(name))]; 
        errorvals = 2*[   std(untreated_in.(name))/sqrt(untreated_in.Ncells-1)...
                        std(peak_in.(name))/sqrt(peak_in.Ncells-1)...
                        std(adapted_in.(name))/sqrt(adapted_in.Ncells-1)]; 
        name = 'ncrfin';
        vals2 = {untreated_in.(name), peak_in.(name), adapted_in.(name)}; 
        means2 = [mean(untreated_in.(name)) mean(peak_in.(name)) mean(adapted_in.(name))]; 
        errorvals2 = 2*[  std(untreated_in.(name))/sqrt(untreated_in.Ncells-1)...
                        std(peak_in.(name))/sqrt(peak_in.Ncells-1)... 
                        std(adapted_in.(name))/sqrt(adapted_in.Ncells-1)]; 
    else
        titlestr = 'LMB';
        yrange = [0 20];

        name = 'rinit';
        vals = {untreatedLMB_in.(name), peakLMB_in.(name), adaptedLMB_in.(name)}; 
        means = [mean(untreatedLMB_in.(name)) mean(peakLMB_in.(name)) mean(adaptedLMB_in.(name))]; 
        errorvals = 2*[   std(untreatedLMB_in.(name))/sqrt(untreatedLMB_in.Ncells-1)...
                        std(peakLMB_in.(name))/sqrt(peakLMB_in.Ncells-1)...
                        std(adaptedLMB_in.(name))/sqrt(adaptedLMB_in.Ncells-1)]; 

        name = 'rfin';
        vals2 = {untreatedLMB_in.(name), peakLMB_in.(name), adaptedLMB_in.(name)}; 
        means2 = [mean(untreatedLMB_in.(name)) mean(peakLMB_in.(name)) mean(adaptedLMB_in.(name))]; 
        errorvals2 = 2*[  std(untreatedLMB_in.(name))/sqrt(untreatedLMB_in.Ncells-1)...
                        std(peakLMB_in.(name))/sqrt(peakLMB_in.Ncells-1)...
                        std(adaptedLMB_in.(name))/sqrt(adaptedLMB_in.Ncells-1)]; 
    end

    clf
    b = bar(cat(1,means,means2)',0.9);%,'FaceColor',[0.1 0.5 0.1]);
    b(1).FaceColor = colors(1,:);
    b(2).FaceColor = colors(2,:);
    hold on
    bins = (1:3);
    errorbar(bins-0.15, means, errorvals,'k','LineStyle','none','linewidth',3);
    errorbar(bins+0.15, means2, errorvals2,'k','LineStyle','none','linewidth',3);
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
    close;
    
    % box plots
    figure, 
    valsvec = [];
    idxvec = [];
    coloridx = [];
    colormap lines;
    lw = 3;
    for i = 1:3
        valsvec = cat(1,valsvec, vals{i}, vals2{i});
        coloridx = cat(1,coloridx, vals{i}*0+1, vals2{i}*0 + 2);
        idxvec = cat(1,idxvec, vals{i}*0+2*i-1, vals2{i}*0+2*i);
    end 
    h = boxplot(valsvec, idxvec,'notch','on','Widths',0.3,...
                            'ColorGroup',coloridx,'Colors',lines(2));
    set(h,'linew',lw);
    set(gca, 'LineWidth', lw);
    set(gcf,'color','w');
    ylim(yrange);
    xticks([1.5 3.5 5.5]);
    set(gca,'XTickLabel', {'ctrl', 'peak', 'adapt'});
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    title(titlestr);
    box off
    axis square
    saveas(gcf, fullfile(dataDir,['rcompare' titlestr '_boxplot.png']));
end

%%
% compare N and C recovery

for LMB = false%[true false]
    
    figure,
    colors = lines(2);
    fs = 30;
    %  LMB = true;

    if ~LMB 
        titlestr = 'no LMB';
        yrange = [];

        name = 'nr';
        vals = {untreated_in.(name) peak_in.(name) adapted_in.(name)};
        means = [mean(untreated_in.(name)) mean(peak_in.(name)) mean(adapted_in.(name))]; 
        errorvals = 2*[   std(untreated_in.(name))/sqrt(untreated_in.Ncells-1)...
                        std(peak_in.(name))/sqrt(peak_in.Ncells-1)...
                        std(adapted_in.(name))/sqrt(adapted_in.Ncells-1)]; 
        name = 'cr';
        vals2 = {untreated_in.(name), peak_in.(name), adapted_in.(name)}; 
        means2 = [mean(untreated_in.(name)) mean(peak_in.(name)) mean(adapted_in.(name))]; 
        errorvals2 = 2*[  std(untreated_in.(name))/sqrt(untreated_in.Ncells-1)...
                        std(peak_in.(name))/sqrt(peak_in.Ncells-1)... 
                        std(adapted_in.(name))/sqrt(adapted_in.Ncells-1)]; 
    else
        titlestr = 'LMB';
        yrange = [];

        name = 'Ninit';
        vals = {untreatedLMB_in.(name), peakLMB_in.(name), adaptedLMB_in.(name)}; 
        means = [mean(untreatedLMB_in.(name)) mean(peakLMB_in.(name)) mean(adaptedLMB_in.(name))]; 
        errorvals = 2*[   std(untreatedLMB_in.(name))/sqrt(untreatedLMB_in.Ncells-1)...
                        std(peakLMB_in.(name))/sqrt(peakLMB_in.Ncells-1)...
                        std(adaptedLMB_in.(name))/sqrt(adaptedLMB_in.Ncells-1)]; 

        name = 'Nfin';
        vals2 = {untreatedLMB_in.(name), peakLMB_in.(name), adaptedLMB_in.(name)}; 
        means2 = [mean(untreatedLMB_in.(name)) mean(peakLMB_in.(name)) mean(adaptedLMB_in.(name))]; 
        errorvals2 = 2*[  std(untreatedLMB_in.(name))/sqrt(untreatedLMB_in.Ncells-1)...
                        std(peakLMB_in.(name))/sqrt(peakLMB_in.Ncells-1)...
                        std(adaptedLMB_in.(name))/sqrt(adaptedLMB_in.Ncells-1)]; 
    end

    clf
    b = bar(cat(1,means,means2)',0.9);%,'FaceColor',[0.1 0.5 0.1]);
    b(1).FaceColor = colors(1,:);
    b(2).FaceColor = colors(2,:);
    hold on
    bins = (1:3);
    errorbar(bins-0.15, means, errorvals,'k','LineStyle','none','linewidth',3);
    errorbar(bins+0.15, means2, errorvals2,'k','LineStyle','none','linewidth',3);
    hold off
    set(gcf,'color','w');
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    box off
    xlim([0.5 max(bins)+0.5]);
    if ~isempty(yrange)
        ylim(yrange);
    end
    axis square
    xticks([1 2 3]);
    set(gca,'XTickLabel', {'ctrl', 'peak', 'adapt'});
    legend({'nr', 'cr'},'Location',legloc{j})
    title(titlestr);
    saveas(gcf, fullfile(dataDir,['Ncompare' titlestr '.png']));
    
    % box plots
    figure, 
    valsvec = [];
    idxvec = [];
    coloridx = [];
    colormap lines;
    lw = 3;
    for i = 1:3
        valsvec = cat(1,valsvec, vals{i}, vals2{i});
        coloridx = cat(1,coloridx, vals{i}*0+1, vals2{i}*0 + 2);
        idxvec = cat(1,idxvec, vals{i}*0+2*i-1, vals2{i}*0+2*i);
    end 
    h = boxplot(valsvec, idxvec,'notch','on','Widths',0.3,...
                            'ColorGroup',coloridx,'Colors',lines(2));
    set(h,'linew',lw);
    set(gca, 'LineWidth', lw);
    set(gcf,'color','w');
    if ~isempty(yrange)
        ylim(yrange);
    end
    xticks([1.5 3.5 5.5]);
    set(gca,'XTickLabel', {'ctrl', 'peak', 'adapt'});
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    title(titlestr);
    box off
    axis square
    saveas(gcf, fullfile(dataDir,['Ncompare' titlestr '_boxplot.png']));
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

Aname = 'A'; % 'Acorr'
instruct = {untreated_in, peak_in, adapted_in};
instructp = {untreatedLMB_in, peakLMB_in, adaptedLMB_in};

for i = 1:3

    input{i} = struct();
        
    input{i}.A = mean(instruct{i}.(Aname));
    input{i}.Ap = mean(instructp{i}.(Aname));
    input{i}.R = mean(instruct{i}.ncrinit);
    input{i}.Rp = mean(instructp{i}.ncrinit);
    input{i}.k = mean(instruct{i}.k);
    input{i}.kp = mean(instructp{i}.k);

    % measured standard errors
    input{i}.sigA = std(instruct{i}.(Aname))/sqrt(instruct{i}.Ncells-1);
    input{i}.sigAp = std(instructp{i}.(Aname))/sqrt(instructp{i}.Ncells-1);
    input{i}.sigR = std(instruct{i}.ncrinit)/sqrt(instruct{i}.Ncells-1);
    input{i}.sigRp = std(instructp{i}.ncrinit)/sqrt(instructp{i}.Ncells-1);
    input{i}.sigk = std(instruct{i}.k)/sqrt(instruct{i}.Ncells-1);
    input{i}.sigkp = std(instructp{i}.k)/sqrt(instructp{i}.Ncells-1);
%     if i == 3
%         input{i}.sigAp = 0.1;
%         input{i}.sigRp = 0.1;
%         input{i}.sigkp = 0.01;
%     end
    
    output{i} = fitKineticModel2(input{i});
end

%% visualize model fit

outputcomb = cat(2,output{:});
inferred = {'kin','kout','cs','ns'};

for i = 1:4
    figure,
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
