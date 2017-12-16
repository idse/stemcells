clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
dataDir = '/Users/idse/data_tmp/0_kinetics/';

FRAPmetadata

%% load FRAP data from different experiments

results = {};
pf = '171214';

for i = 1:numel(FRAPdirs)
	
    fname = fullfile(dataDir,FRAPdirs{i},['results' pf '.mat']);
    
    if exist(fname,'file')
        S = load(fname);
        results{i} = S.allresults;
    end
end

results{1}{1}.good = [0 1];
results{1}{2}.good = [1];

results{2}{3}.good = [1 1];

results{3}{1}.good = [1 1];
results{3}{2}.good = [0 1];
results{3}{3}.good = [1 1];
results{3}{4}.good = [1 1];

results{4}{1}.good = [1 0];

results{6}{1}.good = [1 1];
results{6}{2}.good = [1 1];

results{7}{1}.good = [1 1];
results{7}{2}.good = [1 1 1 1 1 1];

results{8}{2}.good = [1 1];

results{9}{1}.good = [1 1 1];
results{9}{2}.good = [1 1];
results{9}{3}.good = [0 0 1];
results{9}{6}.good = [0 1 1];
results{9}{7}.good = [1 1 0];

results{10}{4}.good = 0*[1 0 1 1 0]; 
results{10}{5}.good = [0 1 1 1]; 

results{12}{1}.good = [0 1];

results{13}{1}.good = [1]; % A50 LMB 2h movie that didn't fully bleach

results{14}{1}.good = [1 1]; % not the best movie

results{15}{1}.good = [1 0];

results{16}{2}.good = [0 1];
results{16}{3}.good = [1 1 1];

results{17}{1}.good = [1 0 0];

results{18}{4}.good = [1 1 1];

results{19}{1}.good = [1 1 1];
results{19}{2}.good = [1 1 1 1 1];
results{19}{3}.good = [0 1 1];
results{19}{4}.good = [1 1];

results{21}{1}.good = [1 1 1 1]; % peak
results{21}{2}.good = [1 1 1 1];
%results{19}{3}.good = [0 0 0 ]; 
results{21}{4}.good = [1 1]; % adapt
results{21}{5}.good = [1 1 1 0 0]; %peak LMB
results{21}{6}.good = [1 1 1 1];
results{21}{8}.good = [1 1 1]; % adapt LMB
results{21}{10}.good = [1 1 1 1 1]; % untreat

for i = 1:10
    results{21}{i}.bc = false;
end

save(fullfile(dataDir, ['FRAPresults' pf '.mat']), 'results');

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

whiskerlength = 1.5;

load(fullfile(dataDir, ['FRAPresults' pf '.mat']));

%short time res: 6 2; 11 4; 
untrIdx = [2 4; 4 1; 6 2; 11 4; 10 5;  18 4; 19 4; 21 10];  % 10 4 is too short
peakIdx = [1 1; 3 1; 3 2; 6 1; 7 1; 12 1; 19 2; 21 1; 21 2]; % 8 2 but some RI mistake
adaptIdx = [3 3; 3 4; 9 2; 7 2; 21 4]; % 19 3 

untrIdxLMB = [2 3; 9 6; 9 7; 11 3; 16 3];
peakIdxLMB = [9 1; 10 2; 11 2; 13 1; 14 1; 19 1; 21 5; 21 6];% 2 1, but that one recoveres very little and too fast
adaptIdxLMB = [9 3; 10 1; 16 2; 17 1; 21 8];

allIdx = {untrIdx, peakIdx, adaptIdx, untrIdxLMB, peakIdxLMB, adaptIdxLMB};
allOut = {};
allLabels = {'untreated', 'peak', 'adapted',...
                    'untreatedLMB', 'peakLMB', 'adaptedLMB'};

debugplots = false;

for i = 1:numel(allIdx)
    legendstr = FRAPdirs(allIdx{i}(:,1));
    allOut{i} = makeInStruct(results, allIdx{i}, debugplots, legendstr, dataDir, allLabels{i});
    close all;
    
    bg = allOut{i}.bg;
%     if i == 5
%         bg = 120;
%     end
    allOut{i}.Araw = (allOut{i}.Nfin - bg)./(allOut{i}.Ninit-bg);
    allOut{i}.ArawCyt = (allOut{i}.Cfin - bg)./(allOut{i}.Cinit-bg);
end

untreated_in = allOut{1};
peak_in = allOut{2}; 
adapted_in = allOut{3}; 

untreatedLMB_in = allOut{4};
peakLMB_in = allOut{5};
adaptedLMB_in = allOut{6};

%% TEST how raw nuclear intensity ratios before and after 
% compare to fitting the recovery amplitude

clf
hold on
disp('---')
for j = 1:6
    
    X = allOut{j};
    bg = 120;%
    %bg = X.bg;
    bla = [X.A X.Araw];

    scatter(bla(:,1),bla(:,2))
    xlabel('fit');
    ylabel('raw');
    axis equal
    %axis([0 0.8 0 0.8]);
    C = corrcoef(bla(:,1),bla(:,2));
    disp(C(1,2))
    %cov(bla(:,1),bla(:,2))
end
legend({'untreat','peak','adapt','untreatLMB','peakLMB','adaptLMB'})
saveas(gcf,fullfile(dataDir,'results', 'fitVraw.png'));

% look at recovery curves to explain lack of correlation in 
% untreated LMB

%%
disp('---')
for j = 5%1:6

    Aall = [];
    Acall = [];
    Nall = [];
    Call = [];

    for i = 8%:size(allIdx{j},1)

        R = results{allIdx{j}(i,1)}{allIdx{j}(i,2)};
        T = R.tracesNuc';

        TC = R.tracesCyt';
        Cinit = mean(TC(1:2,:));
        Cfin = mean(TC(R.tmax - 10:R.tmax,:));
        
        Ninit = mean(T(1:2,:));
        Nfin =  mean(T(R.tmax - 10:R.tmax,:));

        A = R.A(:,1);

        bg = 300.2527;
        Araw = (Nfin - bg)'./(Ninit - bg)';
        Acraw = (Cfin - bg)'./(Cinit - bg)';
        idx = R.good > 0;

        Aall = cat(1,Aall, [A(idx) Araw(idx)]);
        Acall = cat(1,Acall, Acraw);
        
        Nall = cat(1,Nall, [Ninit(idx)', Nfin(idx)']);
        Call = cat(1,Call, [Cinit(idx)', Cfin(idx)']);
    end
    C = corrcoef(Aall(:,1),Aall(:,2));
    disp(C(1,2))
end
%bg = min(min(T(1:10,:)));

%% explain difference

figure,plot(R.tracesNucNorm')
T = R.tracesNuc' - bg;
figure,plot(T./max(T))

%% list all movies in each set 

for i = 1:numel(allIdx)
    
    disp('-----------');
    condIdx = allIdx{i};
    
    for j = 1:size(condIdx,1)
        idx = condIdx(j,:);
        oibfile = oibfiles{idx(1)}{idx(2)};
        disp(sprintf([oibfile '\t\t\t(' FRAPdirs{idx(1)} ')']));
    end
end

%% is measured time scale related to time resolution?

j = 1;
N = size(allIdx{j},1);
x = zeros([1 N]);
y = x; e = x; nt = x; xres = x;
for i = 1:N
    x(i) = mean(results{allIdx{j}(i,1)}{allIdx{j}(i,2)}.k(:,1));
    e(i) = std(results{allIdx{j}(i,1)}{allIdx{j}(i,2)}.k(:,1));
    y(i) = results{allIdx{j}(i,1)}{allIdx{j}(i,2)}.tres;
    nt(i) = size(results{allIdx{j}(i,1)}{allIdx{j}(i,2)}.meanI,1)*y(i);
    xres(i) = results{allIdx{j}(i,1)}{allIdx{j}(i,2)}.xyres;
end
u = nt;
[~,sorti] = sort(u);
errorbar(u(sorti),x(sorti),e(sorti),'-x')

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
    errorvals = 4*[   std(untreated_in.(name))/sqrt(untreated_in.Ncells-1)...
                    std(peak_in.(name))/sqrt(peak_in.Ncells-1)...
                    std(adapted_in.(name))/sqrt(adapted_in.Ncells-1)]; 
                %[std(untreated_in.(name)) std(peak_in.(name)) std(adapted_in.(name))]; 
    
    vals2 = {untreatedLMB_in.(name), peakLMB_in.(name), adaptedLMB_in.(name)}; 
    means2 = [mean(untreatedLMB_in.(name)) mean(peakLMB_in.(name)) mean(adaptedLMB_in.(name))]; 
    errorvals2 = 4*[  std(untreatedLMB_in.(name))/sqrt(untreatedLMB_in.Ncells-1)...
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
    h = boxplot(valsvec, idxvec,'notch','off','Widths',0.3,...
                            'ColorGroup',coloridx,'Colors',lines(2),...
                            'Whisker',whiskerlength); % 1.5 default
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
                            'ColorGroup',coloridx,'Colors',lines(2),...
                            'Whisker',whiskerlength); 
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

for LMB = [true false]%false
    
    figure,
    colors = lines(2);
    fs = 30;
    %  LMB = true;

    if ~LMB 
        titlestr = 'recovery';
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
        titlestr = 'recovery LMB';
        yrange = [];

        name = 'nr';
        vals = {untreatedLMB_in.(name), peakLMB_in.(name), adaptedLMB_in.(name)}; 
        means = [mean(untreatedLMB_in.(name)) mean(peakLMB_in.(name)) mean(adaptedLMB_in.(name))]; 
        errorvals = 2*[   std(untreatedLMB_in.(name))/sqrt(untreatedLMB_in.Ncells-1)...
                        std(peakLMB_in.(name))/sqrt(peakLMB_in.Ncells-1)...
                        std(adaptedLMB_in.(name))/sqrt(adaptedLMB_in.Ncells-1)]; 

        name = 'cr';
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
    ylim([0 1]);
    xticks([1 2 3]);
    set(gca,'XTickLabel', {'ctrl', 'peak', 'adapt'});
    legend({'nuclear', 'cytoplasmic'},'Location','NorthEast','FontSize',20)
    title(titlestr);
    axis square
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
                            'ColorGroup',coloridx,'Colors',lines(2),...
                            'Whisker',whiskerlength);
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

Aname = 'Araw'; % 'Acorr'
instruct = {untreated_in, peak_in, adapted_in};
instructp = {untreatedLMB_in, peakLMB_in, adaptedLMB_in};

disp('-------------------------------------------------');

for i = 1:3

    input{i} = struct();

    input{i}.A = mean(instruct{i}.(Aname));
    input{i}.Ap = mean(instructp{i}.(Aname));
    input{i}.R = mean(instruct{i}.ncrinit);
    input{i}.Rp = mean(instructp{i}.ncrinit);
    input{i}.k = mean(instruct{i}.k);
    input{i}.kp = mean(instructp{i}.k);
    
    input{i}.Ac = mean(instruct{i}.ArawCyt);
    input{i}.Acp = mean(instructp{i}.ArawCyt);
    
    % measured standard errors
    input{i}.sigA = std(instruct{i}.(Aname))/sqrt(instruct{i}.Ncells-1);
    input{i}.sigAp = std(instructp{i}.(Aname))/sqrt(instructp{i}.Ncells-1);
    input{i}.sigR = std(instruct{i}.ncrinit)/sqrt(instruct{i}.Ncells-1);
    input{i}.sigRp = std(instructp{i}.ncrinit)/sqrt(instructp{i}.Ncells-1);
    input{i}.sigk = std(instruct{i}.k)/sqrt(instruct{i}.Ncells-1);
    input{i}.sigkp = std(instructp{i}.k)/sqrt(instructp{i}.Ncells-1);
    
    input{i}.sigAc = std(instruct{i}.ArawCyt)/sqrt(instruct{i}.Ncells-1);
    input{i}.sigAcp = std(instructp{i}.ArawCyt)/sqrt(instructp{i}.Ncells-1);
end

alpha = 1.13;
output = fitKineticModel4(input);%,alpha);
outputNew = fitKineticModelNew(input);%,alpha);

%%
output2 = fitKineticModel3(input,alpha);

%%
output2 = fitKineticModelFixKin(input);

%%
output2 = fitKineticModelFixKinNoCytSeq(input);

%% consistency check 

% changing R
output = output2;
alpha = output2{1}.alpha;

kap = @(kin, kout) kin/(kin+kout);
A = @(kin, kout, cs, ns) kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs)/(ns + kap(kin,kout)*(1-ns-cs));
R = @(kin, kout, cs, ns) alpha*(kout*ns + kin*(1-cs))/(kin*cs + kout*(1-ns));
k = @(kin, kout) kin + kout;

A = @(kin, kout, cs, ns) kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs)/(ns + kap(kin,kout)*(1-ns-cs));
R = @(kin, kout, cs, ns) alpha*(ns + kap(kin,kout)*(1-ns-cs))/(cs + (1-kap(kin,kout))*(1-ns-cs));

for i = 1:3
    N = 2;
    disp('--------------------');
    disp(['A:  ' num2str([A(output{i}.kin, output{i}.kout, output{i}.cs, output{i}.ns) input{i}.A],N)]);
    disp(['R:  ' num2str([R(output{i}.kin, output{i}.kout, output{i}.cs, output{i}.ns) input{i}.R],N)]);
    disp(['k:  ' num2str([k(output{i}.kin, output{i}.kout) input{i}.k],N)]);
    disp(['Ap: ' num2str([A(output{i}.kin, output{i}.koutp, output{i}.cs, output{i}.ns) input{i}.Ap],N)]);
    disp(['Rp: ' num2str([R(output{i}.kin, output{i}.koutp, output{i}.cs, output{i}.ns) input{i}.Rp],N)]);
    disp(['kp: ' num2str([k(output{i}.kin, output{i}.koutp) input{i}.kp],N)]);
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
    errorvals = 4*[outputcomb.(errname)];
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

%% visualize model fit differently

outputcomb = cat(2,output{:});
inferred = {'kin','kout','cs','ns'};
vals = {};
errorvals = {};

figure,
hold on
%s = 10^3;
s=1;
k=1;
for i = 3:4
    name = inferred{i};
    errname = ['sig' name];
    vals{k} = s*2*[outputcomb.(name)];
    errorvals{k} = s*2*[outputcomb.(errname)];
    k = k+1;
end
hold off
   
fs = 26;
b = bar(cat(1,vals{2}, vals{1})',0.9);%'FaceColor',[0.1 0.5 0.1]);
colors = lines(2);
b(1).FaceColor = colors(1,:);
b(2).FaceColor = colors(2,:);
hold on
bins = (1:3);
errorbar(bins - 0.15, vals{2}, errorvals{2},...
                        'k','linestyle','none','linewidth',2);
errorbar(bins + 0.15, vals{1}, errorvals{1},...
                        'k','linestyle','none','linewidth',2);
hold off
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
box off
set(gca,'XTickLabel', {'untreated', 'peak', 'adapt'});

%ylabel('10^{-3} sec^{-1}');
%legend({'k_{out}','k_{in}'});
%title('nuclear exchange')
%fname = 'nuclearExchange_values.png';

ylabel('10^{-3} sec^{-1}');
legend({'nuclear','cytoplasm'});
title('sequestration')
fname = 'sequestered_values.png';
ylim([0 0.5]);

xlim([0.5 3.5]);
axis square
saveas(gcf, fullfile(dataDir, fname));

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
