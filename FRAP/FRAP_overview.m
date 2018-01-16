clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
dataDir = '/Users/idse/data_tmp/0_kinetics/';

FRAPmetadata

%% load FRAP data from different experiments

results = {};
pf = '171217c';%'171214';

for i = 1:numel(FRAPdirs)

    fname = fullfile(dataDir,FRAPdirs{i},['results' pf '.mat']);

    if exist(fname,'file')
        S = load(fname);
        results{i} = S.allresults;
        for j = 1:numel(results{i})
            results{i}{j}.frapframe = frapframes(i);
            if i==11 && j==2
                % an exception to the rule that frapframe is the same for
                % all movies in a set
                results{i}{j}.frapframe = 5;
            end
        end
    end
end

results{1}{1}.good = [0 1];
results{1}{2}.good = [1];

results{2}{3}.good = [1 1];

results{3}{1}.good = [1 1];
results{3}{2}.good = [1 1];
results{3}{3}.good = [0 1];
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

results{11}{4}.good = [1 1 0 1 1];

results{12}{1}.good = [0 1];

results{13}{1}.good = [1]; % A50 LMB 2h movie that didn't fully bleach

results{14}{1}.good = [1 0]; % not the best movie

results{15}{1}.good = [1 0];

results{16}{2}.good = [0 1];
results{16}{3}.good = [1 1 1];

results{17}{1}.good = [1 0 0];

results{18}{4}.good = [0 1 0];

results{19}{1}.good = [1 1 1];
results{19}{2}.good = [1 1 1 1 1];
results{19}{3}.good = [1 0 1];
results{19}{4}.good = [1 1];

results{21}{1}.good = [1 1 1 1]; % peak
results{21}{2}.good = [1 1 1 1];
%results{19}{3}.good = [0 0 0 ]; 
results{21}{4}.good = [1 1]; % adapt
results{21}{5}.good = [1 0 1 0 1]; %peak LMB
results{21}{6}.good = [1 1 1 1];
results{21}{8}.good = [1 1 1]; % adapt LMB
results{21}{10}.good = [0 1 1 1 1]; % untreat CHANGED

results{22}{3}.good = [0 1]; % adapted
results{22}{4}.good = [1 1 1]; % adapted

results{23}{5}.good = [1 1 1 0]; % LMB adapted
results{23}{6}.good = [1 1 1]; % LMB adapted

save(fullfile(dataDir, ['FRAPresults' pf '.mat']), 'results');


% for i = 1:numel(FRAPdirs)
% 	
%     fname = fullfile(dataDir,FRAPdirs{i},['results' pf '.mat']);
%     if exist(fname,'file')
%         allresults = results{i};
%         save(fname, 'allresults');
%     end
% end

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

load(fullfile(dataDir, ['FRAPresults' pf '.mat']));

%short time res: 6 2; 11 4; NOW 10 5; 
untrIdx = [6 2; 11 4; 18 4; 19 4; 21 10]; % 10 5 good or bad
%4 1; 
% 10 4 is too short
% 2 4 is outlier in terms of speed, but also taken with totally different
% imaging conditions (zoom), so excluded bc perhaps we're seeing strong diffusive
% recovery
peakIdx = [1 1; 3 2; 6 1; 7 1; 19 2; 21 1]; 
% 21 2
% 8 2 but some RI mistake
% 12 1; looks pretty adapted already
adaptIdx = [3 4; 9 2; 19 3; 21 4; 22 3; 22 4]; % 19 3  
% 7 2 starts with n:c too high, seems not fully adapted, and recovers less
% 3 3 drifts out of focus, ok for nucleus but makes cytoplasmic 
% readout bad right away

%untrIdxLMB = [2 3; 9 6; 9 7; 11 3; 16 3];
%peakIdxLMB = [9 1; 10 2; 11 2; 13 1; 14 1; 19 1; 21 5; 21 6];% 2 1, but that one recoveres very little and too fast
untrIdxLMB = [2 3; 9 6; 9 7; 11 2; 11 3; 16 3];
peakIdxLMB = [9 1; 10 2; 11 1; 14 1; 21 5;];% 2 1, but that one recoveres very little and too fast
% 13 1; garbage trace: too much movement?
% 19 1; LMB not long enough, recovery too fast, bad cyto
%  21 6
adaptIdxLMB = [9 3; 10 1; 17 1; 21 8; 23 5; 23 6]; %16 2; 

allIdx = {untrIdx, peakIdx, adaptIdx, untrIdxLMB, peakIdxLMB, adaptIdxLMB};
allOut = {};
allLabels = {'untreated', 'peak', 'adapted',...
                    'untreatedLMB', 'peakLMB', 'adaptedLMB'};

debugplots = false;

for i = 1:numel(allIdx)
    legendstr = FRAPdirs(allIdx{i}(:,1));
    allOut{i} = makeInStruct(results, allIdx{i}, debugplots, legendstr, dataDir, allLabels{i},bleachCorrect);
    
    bg = allOut{i}.bg;
%     if i == 5
%         bg = 120;
%     end
    allOut{i}.Araw = (allOut{i}.Nfin - allOut{i}.Nbleach)./(allOut{i}.Ninit-bg);
    allOut{i}.Braw = (allOut{i}.Nbleach - bg)./(allOut{i}.Ninit-bg);
    allOut{i}.bn = allOut{i}.Braw;
    
    allOut{i}.ArawCyt = (allOut{i}.Cbleach-allOut{i}.Cfin)./(allOut{i}.Cinit-bg);
    allOut{i}.BrawCyt = (allOut{i}.Cfin - bg)./(allOut{i}.Cinit-bg);
    allOut{i}.bc = (allOut{i}.Cbleach - bg)./(allOut{i}.Cinit-bg);
    
    allOut{i}.Acorr = allOut{i}.Araw./(allOut{i}.bc - allOut{i}.bn);
    allOut{i}.AcorrCyt = allOut{i}.ArawCyt./(allOut{i}.bc - allOut{i}.bn);
    
    %allOut{i}.Rfin = (allOut{i}.Nfin-bg)./(allOut{i}.Cfin-bg);
    %allOut{i}.Rinit = (allOut{i}.Ninit-bg)./(allOut{i}.Cinit-bg);
    allOut{i}.Rrat = allOut{i}.ncrfin./allOut{i}.ncrinit;
    
    allOut{i}.tau = 1./(60*allOut{i}.k);
end

untreated_in = allOut{1};
peak_in = allOut{2}; 
adapted_in = allOut{3}; 

untreatedLMB_in = allOut{4};
peakLMB_in = allOut{5};
adaptedLMB_in = allOut{6};

%% display relevant values for cleanup

names = {'Araw','Braw','ArawCyt','BrawCyt','bc','Acorr','AcorrCyt','k',...
            'Rinit','Rfin','Rrat','tau'};
fs = 15;
for i = 8%[1 3 8]%9:10%1:8
    name = names{i};
    disp(['--' name '--']);
    figure,
    hold on
    for j = 1:3
        N = numel(allOut{j}.(name));
        scatter(allOut{j}.(name)*0 + j + 0.1, allOut{j}.(name),'LineWidth',2)
        disp([nanmean(allOut{j}.(name)) nanstd(allOut{j}.(name)) 2*nanstd(allOut{j}.(name))/sqrt(N-1)])
        for k = 1:N
            text(j+0.1, allOut{j}.(name)(k),...
                [FRAPdirs{allOut{j}.movieIdx(k,1)}(3:6) '.' num2str(allOut{j}.movieIdx(k,:))],'LineWidth',2)
        end
    end
    disp('--');
    for j = 1:3
        N = numel(allOut{j+3}.k);
        scatter(allOut{j+3}.(name)*0 + j+1/2, allOut{j+3}.(name),'LineWidth',2)
        disp([mean(allOut{j+3}.(name)) std(allOut{j+3}.(name)) 2*std(allOut{j+3}.(name))/sqrt(N-1)])
        for k = 1:N
            text(j+1/2+0.05, allOut{j+3}.(name)(k),...
                [FRAPdirs{allOut{j+3}.movieIdx(k,1)}(3:6) '.' num2str(allOut{j+3}.movieIdx(k,:))],'LineWidth',2)
        end
    end
    hold off
    title(name);
    set(gcf,'color','w');
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    saveas(gcf,fullfile(dataDir, ['measuredVals_' names{i} '.png']));
    %close;
end

%% TEST how raw nuclear intensity ratios before and after 
% compare to fitting the recovery amplitude

% minor inconsistencies come from the fact that readFRAPprofile2 uses
% frapframe = 3 instead of the actual frapframe for determining the
% normalization

figure,
clf
hold on
disp('---')

%untreated [18 4] doesn't correlate because normalized to max and max is
%not prebleach but end because of drift there:
% Araw is correct

for j = 1:6
    
    X = allOut{j};
    bla = [X.A X.Araw];

    scatter(bla(:,1),bla(:,2))
    B = [bla(:,1),bla(:,2)]
    xlabel('fit');
    ylabel('raw');
    axis equal
    axis([0 0.8 0 0.8]);
    C = corrcoef(bla(:,1),bla(:,2));
    disp(C(1,2))
    %disp([num2str([mean(X.Araw) std(X.Araw)]) '     fit: ' num2str([mean(X.A) std(X.A)])])
    %cov(bla(:,1),bla(:,2))
end
legend({'untreat','peak','adapt','untreatLMB','peakLMB','adaptLMB'})
saveas(gcf,fullfile(dataDir,'results', 'fitVrawB.png'));

% look at recovery curves to explain lack of correlation in 
% untreated LMB

% exponential in fit starts at 0, not at FRAP frame?~!

%% cytoplasmic bleach analysis

cytresults = {};
pf = '';

for i = 1:numel(FRAPdirs)

    fname = fullfile(dataDir,FRAPdirs{i},['resultsCyt' pf '.mat']);

    if exist(fname,'file')
        S = load(fname);
        cytresults{i} = S.allresults;
        for j = 1:numel(cytresults{i})
            cytresults{i}{j}.frapframe = frapframes(i);
            if i==26 && j==3
                % an exception to the rule that frapframe is the same for
                % all movies in a set
                cytresults{i}{j}.frapframe = 2;
            end
        end
    end
end

cytresults{27}{2}.good = [0 1];
cytresults{25}{1}.good = [0 1 1];

% combine similar condition measurements for cytoplasm

% 26 6; = 0612 untreated3, weird curves, may be ok after redef masks
% 29 1; 30 1 Untr_Nucl-Cyt_Slow folders, much cell deformation
% 26 4; is ok but outlier
untrIdx = [12 1; 25 4; 26 5; 27 3; 27 4; 28 1];
peakIdx = [25 3; 25 1; 25 2; 26 1; 27 1];
adaptIdx = [26 2; 26 3; 27 2];

untrIdxLMB = [14 1];
peakIdxLMB = [];
adaptIdxLMB = [];

allIdxCyt = {untrIdx, peakIdx, adaptIdx, untrIdxLMB, peakIdxLMB, adaptIdxLMB};

allLabels = {'untreated', 'peak', 'adapted',...
                    'untreatedLMB', 'peakLMB', 'adaptedLMB'};
debugplots = false;

for i = 1:3
    legendstr = FRAPdirs(allIdxCyt{i}(:,1));
    
    allOutCyt{i} = makeInStruct(cytresults, allIdxCyt{i}, debugplots, legendstr, dataDir, allLabels{i},bleachCorrect);
    
    bg = allOutCyt{i}.bg;
    
    allOutCyt{i}.Araw = -(allOutCyt{i}.Nfin - allOutCyt{i}.Nbleach)./(allOutCyt{i}.Ninit-bg);
    allOutCyt{i}.Braw = (allOutCyt{i}.Nfin - bg)./(allOutCyt{i}.Ninit-bg);
    allOutCyt{i}.bn = (allOutCyt{i}.Nbleach - bg)./(allOutCyt{i}.Ninit-bg);
    
    allOutCyt{i}.ArawCyt = (allOutCyt{i}.Cfin - allOutCyt{i}.Cbleach)./(allOutCyt{i}.Cinit-bg);
    allOutCyt{i}.BrawCyt = (allOutCyt{i}.Cbleach - bg)./(allOutCyt{i}.Cinit-bg);
	allOutCyt{i}.bc = (allOutCyt{i}.Cbleach - bg)./(allOutCyt{i}.Cinit-bg);
    
    allOutCyt{i}.Acorr = allOutCyt{i}.Araw./(allOutCyt{i}.bn - allOutCyt{i}.bc);
    
    allOutCyt{i}.Rfin = (allOutCyt{i}.Nfin-bg)./(allOutCyt{i}.Cfin-bg);
    allOutCyt{i}.Rinit = (allOutCyt{i}.Ninit-bg)./(allOutCyt{i}.Cinit-bg);
    allOutCyt{i}.Rrat = allOutCyt{i}.Rfin./allOutCyt{i}.Rinit;
end

%
% 
% %%
% i = 1
% r = (allOutCyt{i}.Cfin - allOutCyt{i}.Cbleach)./(allOutCyt{i}.Cinit-allOutCyt{i}.bg);
% [r allOutCyt{i}.Cfin allOutCyt{i}.Cbleach allOutCyt{i}.Cinit allOutCyt{i}.bg]

%%
disp('compare nuclear:cytoplasmic recovery');

titlestr = 'rcompare';

figure,
hold on
    
vals = {};
vals2 = {};

for i = 1:3
    
    bg = allOutCyt{i}.bg;
    %allOutCyt{i}.ncrfin./allOutCyt{i}.ncrinit;
    [allOutCyt{i}.Nfin allOutCyt{i}.Cfin bg];
    
    Rfin = (allOutCyt{i}.Nfin-bg)./(allOutCyt{i}.Cfin-bg);
    Rinit = (allOutCyt{i}.Ninit-bg)./(allOutCyt{i}.Cinit-bg);
    %[Rinit Rfin Rfin./Rinit]
    Rinit_CBmean = mean(Rinit);
    Rfin_CBmean = mean(Rfin);
    Rrat_CBmean = mean(Rfin./Rinit);
    
    vals{i} = allOut{i}.Rrat;
    vals2{i} = allOutCyt{i}.Rrat;
    
    %[allOut{i}.ncrinit allOut{i}.ncrfin allOut{i}.ncrfin./allOut{i}.ncrinit]
    Rinit_NBmean = nanmean(allOut{i}.ncrinit);
    Rfin_NBmean = nanmean(allOut{i}.ncrfin);
    Rrat_NBmean = nanmean(allOut{i}.Rrat);
    
    scatter(allOut{i}.Rrat*0 + i, allOut{i}.Rrat);
    scatter(allOutCyt{i}.Rrat*0 + i+0.5, allOutCyt{i}.Rrat);
    
    disp(num2str([Rinit_CBmean Rinit_NBmean Rfin_CBmean Rfin_NBmean Rrat_CBmean Rrat_NBmean],2))
end

%% display relevant values for cleanup cytoplasm

% k is nuclear recovery time??
% A is the nuclear amplitude

names = {'Araw','Braw','ArawCyt','BrawCyt','bc','Acorr','AcorrCyt','k',...
            'Rinit','Rfin','Rrat','tau'};
fs = 15;
for i = 11%[1 3 8]%9:10%1:8
    name = names{i};
    disp(['--' name '--']);
    figure,
    hold on
    for j = 1:3
        N = numel(allOutCyt{j}.(name));
        scatter(allOutCyt{j}.(name)*0 + j + 0.1, allOutCyt{j}.(name),'LineWidth',2)
        disp([nanmean(allOutCyt{j}.(name)) nanstd(allOutCyt{j}.(name)) 2*nanstd(allOutCyt{j}.(name))/sqrt(N-1)])
        for k = 1:N
            text(j+0.1, allOutCyt{j}.(name)(k),...
                [FRAPdirs{allOutCyt{j}.movieIdx(k,1)}(3:6) '.' num2str(allOutCyt{j}.movieIdx(k,:))],'LineWidth',2)
        end
    end
    hold off
    title(name);
    set(gcf,'color','w');
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    saveas(gcf,fullfile(dataDir, ['measuredVals_' names{i} '.png']));
    %close;
end

%% box plots for R'/R

figure, 
hold on
valsvec = [];
idxvec = [];
coloridx = [];
colormap lines;
lw = 3;
fs = 25;
whiskerlength = Inf;
for i = 1:3
    valsvec = cat(1,valsvec, vals{i}, vals2{i});
    coloridx = cat(1,coloridx, vals{i}*0+1, vals2{i}*0 + 2);
    idxvec = cat(1,idxvec, vals{i}*0+2*i-1, vals2{i}*0+2*i);
end 
h = boxplot(valsvec, idxvec,'notch','off','Widths',0.5,...
                        'ColorGroup',coloridx,'Colors',lines(2),...
                        'Whisker',whiskerlength); 
set(h,'linew',lw);
set(gca, 'LineWidth', lw);
set(gcf,'color','w');
ylim([0 3]);
xticks([1.5 3.5 5.5]);
set(gca,'XTickLabel', {'untreated', 'peak', 'adapt'});
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
%title('change');
box off
plot([0 8], [1 1],'k--','LineWidth',2);
pbaspect([3 2 1]);
%axis square
hold off
saveas(gcf, fullfile(dataDir,['rcompare' titlestr '_boxplot.png']));

%% box plots for R'/R just nuclear

figure, 
hold on
valsvec = [];
idxvec = [];
coloridx = [];
colormap lines;
lw = 3;
fs = 25;
whiskerlength = Inf;
for i = 1:3
    valsvec = cat(1,valsvec, vals{i});
    coloridx = cat(1,coloridx, vals{i}*0+1);
    idxvec = cat(1,idxvec, vals{i}*0+2*i-1);
end 
h = boxplot(valsvec, idxvec,'notch','off','Widths',0.5,...
                        'ColorGroup',coloridx,'Colors',lines(2),...
                        'Whisker',whiskerlength); 
set(h,'linew',lw);
set(gca, 'LineWidth', lw);
set(gcf,'color','w');
ylim([0 1]);
xticks([1 2 3]);
set(gca,'XTickLabel', {'untreated', 'peak', 'adapt'});
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
%title('change');
box off
pbaspect([3 2 1]);
%axis square
hold off
saveas(gcf, fullfile(dataDir,['rcompareJustNuclear' titlestr '_boxplot.png']));


%%
%LMB
for i = 4:6
    
    %[allOut{i}.ncrinit allOut{i}.ncrfin allOut{i}.ncrfin./allOut{i}.ncrinit]
    Rinit_NBmean = nanmean(allOut{i}.ncrinit);
    Rfin_NBmean = nanmean(allOut{i}.ncrfin);
    Rrat_NBmean = nanmean(allOut{i}.Rrat);
    
    scatter(allOut{i}.Rrat*0 + i, allOut{i}.Rrat);
    
    disp(num2str([ Rinit_NBmean  Rfin_NBmean  Rrat_NBmean],2))
end

%[allOutCyt{i}.Nfin allOutCyt{i}.Cfin bg Rfin];

%% list all movies in each set 

for i = 1:numel(allIdx)
    
    disp('------------');
    condIdx = allIdx{i};
    
    for j = 1:size(condIdx,1)
        idx = condIdx(j,:);
        oibfile = oibfiles{idx(1)}{idx(2)};
        fprintf([num2str(idx(1)) ' ' num2str(idx(2)) '\t' oibfile '(' num2str(sum(results{idx(1)}{idx(2)}.good)) ') \t\t\t(' FRAPdirs{idx(1)} ')\n']);
    end
end

%% display amounts of input data

Nmovies = [untreated_in.Nmovies peak_in.Nmovies adapted_in.Nmovies]
NmoviesLMB = [untreatedLMB_in.Nmovies peakLMB_in.Nmovies adaptedLMB_in.Nmovies]

Ncells = [sum(untreated_in.Ncells) sum(peak_in.Ncells) sum(adapted_in.Ncells)]
NcellsLMB = [sum(untreatedLMB_in.Ncells) sum(peakLMB_in.Ncells) sum(adaptedLMB_in.Ncells)]


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


%%
whiskerlength = Inf;

disp('display measured parameters');

measured = {'A','Acorr','k','x'};% 'rinit','rfin', % w and w/o LMB
legloc = {'NorthWest','NorthEast','NorthEast','NorthWest'};
yranges = {[0 0.8],[0 1],[0 0.008]};%[0 20],[0 5]};
titles = {'recovery amplitude A','','recovery rate k'};
for j = 1:3%1:3

    name = measured{j};
    
    vals = {untreated_in.(name), peak_in.(name), adapted_in.(name)}; 
    means = [mean(untreated_in.(name)) mean(peak_in.(name)) mean(adapted_in.(name))]; 
    errorvals = 2*[   std(untreated_in.(name))/sqrt(sum(untreated_in.Ncells)-1)...
                    std(peak_in.(name))/sqrt(sum(peak_in.Ncells)-1)...
                    std(adapted_in.(name))/sqrt(sum(adapted_in.Ncells)-1)]; 
                %[std(untreated_in.(name)) std(peak_in.(name)) std(adapted_in.(name))]; 
    
    vals2 = {untreatedLMB_in.(name), peakLMB_in.(name), adaptedLMB_in.(name)}; 
    means2 = [mean(untreatedLMB_in.(name)) mean(peakLMB_in.(name)) mean(adaptedLMB_in.(name))]; 
    errorvals2 = 2*[  std(untreatedLMB_in.(name))/sqrt(sum(untreatedLMB_in.Ncells)-1)...
                    std(peakLMB_in.(name))/sqrt(sum(peakLMB_in.Ncells)-1)...
                    std(adaptedLMB_in.(name))/sqrt(sum(adaptedLMB_in.Ncells)-1)]; 
                %[std(untreatedLMB_in.(name)) std(peakLMB_in.(name)) std(adaptedLMB_in.(name))]; 

%     figure, 
%     colors = lines(2);
%     fs = 30;
%     b = bar(cat(1,means,means2)',0.9);%,'FaceColor',[0.1 0.5 0.1]);
%     b(1).FaceColor = colors(1,:);
%     b(2).FaceColor = colors(2,:);
%     hold on
%     errorbar([1 2 3]-0.15,means, errorvals,'k','LineStyle','none','linewidth',3);
%     errorbar([1 2 3]+0.15,means2, errorvals2,'k','LineStyle','none','linewidth',3);
%     hold off
%     set(gcf,'color','w');
%     set(gca, 'LineWidth', 2);
%     set(gca,'FontSize', fs)
%     set(gca,'FontWeight', 'bold')
%     box off
%     xlim([0.5 3.5]);
%     if ~isempty(yranges{j})
%         ylim(yranges{j});
%     end
%     axis square
%     xticks([1 2 3]);
%     set(gca,'XTickLabel', {'ctrl', 'peak', 'adapt'});
%     legend({name, [name '+LMB']},'Location',legloc{j})
    %saveas(gcf, fullfile(dataDir,[name '_bargraph.png']));
    
    
%     % box plots
%     figure, 
%     valsvec = [];
%     idxvec = [];
%     coloridx = [];
%     colormap lines;
%     lw = 3;
%     for i = 1:3
%         valsvec = cat(1,valsvec, vals{i}, vals2{i});
%         coloridx = cat(1,coloridx, vals{i}*0+1, vals2{i}*0 + 2);
%         idxvec = cat(1,idxvec, vals{i}*0+2*i-1, vals2{i}*0+2*i);
%     end
%     h = boxplot(valsvec, idxvec,'notch','on','Widths',0.3,...
%                             'ColorGroup',coloridx,'Colors',lines(2),...
%                             'Whisker',whiskerlength); % 1.5 default
%     set(h,'linew',lw);
%     set(gca, 'LineWidth', lw);
%     set(gcf,'color','w');
%     if ~isempty(yranges{j})
%         ylim(yranges{j});
%     end
%     xticks([1.5 3.5 5.5]);
%     set(gca,'XTickLabel', {'ctrl', 'peak', 'adapt'});
%     set(gca,'FontSize', fs)
%     set(gca,'FontWeight', 'bold')
%     box off
%     axis square
%     saveas(gcf, fullfile(dataDir,[name '_boxplot.png']));

    % reduced box plots
    figure, 
    valsvec = [];
    idxvec = [];
    coloridx = [];
    colormap lines;
    lw = 3;
    for i = 1:3
        valsvec = cat(1,valsvec, vals{i});
        coloridx = cat(1,coloridx, i);%vals{i}*0);
        idxvec = cat(1,idxvec, vals{i}*0+2*i-1);
    end
    h = boxplot(valsvec, idxvec,'notch','off','Widths',0.5,...
                            'ColorGroup',coloridx,'Colors',lines(3),...
                            'Whisker',whiskerlength); % 1.5 default
    set(h,'linew',lw);
    set(gca, 'LineWidth', lw);
    set(gcf,'color','w');
    if ~isempty(yranges{j})
        ylim(yranges{j});
    end
    %title(titles{j});
    xticks([1 2 3]);
    set(gca,'XTickLabel', {'untreat', 'peak', 'adapt'});
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    box off
    axis square
    saveas(gcf, fullfile(dataDir,[name '_boxplotreduced.png']));
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
        errorvals = 2*[   std(untreated_in.(name))/sqrt(sum(untreated_in.Ncells)-1)...
                    std(peak_in.(name))/sqrt(sum(peak_in.Ncells)-1)...
                    std(adapted_in.(name))/sqrt(sum(adapted_in.Ncells)-1)]; 
                
        name = 'ncrfin';
        vals2 = {untreated_in.(name), peak_in.(name), adapted_in.(name)}; 
        means2 = [mean(untreated_in.(name)) mean(peak_in.(name)) mean(adapted_in.(name))]; 
        errorvals2 = 2*[   std(untreated_in.(name))/sqrt(sum(untreated_in.Ncells)-1)...
                    std(peak_in.(name))/sqrt(sum(peak_in.Ncells)-1)...
                    std(adapted_in.(name))/sqrt(sum(adapted_in.Ncells)-1)]; 
    else
        titlestr = 'LMB';
        yrange = [0 20];

        name = 'rinit';
        vals = {untreatedLMB_in.(name), peakLMB_in.(name), adaptedLMB_in.(name)}; 
        means = [mean(untreatedLMB_in.(name)) mean(peakLMB_in.(name)) mean(adaptedLMB_in.(name))]; 
        errorvals = 2*[  std(untreatedLMB_in.(name))/sqrt(sum(untreatedLMB_in.Ncells)-1)...
                    std(peakLMB_in.(name))/sqrt(sum(peakLMB_in.Ncells)-1)...
                    std(adaptedLMB_in.(name))/sqrt(sum(adaptedLMB_in.Ncells)-1)]; 

        name = 'rfin';
        vals2 = {untreatedLMB_in.(name), peakLMB_in.(name), adaptedLMB_in.(name)}; 
        means2 = [mean(untreatedLMB_in.(name)) mean(peakLMB_in.(name)) mean(adaptedLMB_in.(name))]; 
        errorvals2 = 2*[  std(untreatedLMB_in.(name))/sqrt(sum(untreatedLMB_in.Ncells)-1)...
                    std(peakLMB_in.(name))/sqrt(sum(peakLMB_in.Ncells)-1)...
                    std(adaptedLMB_in.(name))/sqrt(sum(adaptedLMB_in.Ncells)-1)]; 
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

%untrIdx = [6 2; 11 4; 18 4; 19 4; 21 10];  
%peakIdx = [1 1; 3 2; 6 1; 7 1; 19 2; 21 1]; 
%adaptIdx = [3 4; 9 2; 19 3; 21 4; 22 3; 22 4]; % 19 3  

colors = lines(3);
idx = {[2 4], [1 1], [3 4]};
shapeIndices = [1 1 2];

shift = [5 0 -2];

% NEXT: change fit here to get nice fit for figure
clf
hold on
h = [];
for i = 1:3

    ei = idx{i}(1);
    fi = idx{i}(2);
    shapeIdx = shapeIndices(i);
    disp([FRAPdirs{ei} '/' results{ei}{fi}.description])

    tracesnorm = results{ei}{fi}.tracesNucNorm;%(shapeIdx,:);
    t = results{ei}{fi}.tres*(0:size(tracesnorm,2)-1);
    tmax = results{ei}{fi}.tmax;
    tlim = round(max(tmax)*results{ei}{fi}.tres);
    frapframe = results{ei}{fi}.frapframe;

    t = t + shift(i)*results{ei}{fi}.tres;
    lw = 1;
    h(i) = plot(t, tracesnorm(shapeIdx,:),...
                            'Color',colors(i,:),'LineWidth',lw)

    fitval = fitFRAP2(results{ei}{fi}); % change to 3 for recent version
    A = fitval.A(shapeIdx,1);
    k = fitval.k(shapeIdx,1);
    B = fitval.B(shapeIdx,1);

    func = @(p1,p2,p3,x) p1*(1-exp(-x*p2)) + p3;
    fitcurve = func(A,k,B,t(frapframe:end)-t(frapframe));  
    
    lw = 3;
    plot(t(frapframe:end),fitcurve,...
                            'Color', colors(i,:),'LineWidth',lw)
end

hold off

fs = 26;
%ylabel('relative recovery', 'FontSize',fs, 'FontWeight','Bold');
%xlabel('time (sec)', 'FontSize',fs, 'FontWeight','Bold');
xlim([0 1000]);
ylim([0 0.6]);
set(gcf,'color','w');
set(gca, 'LineWidth', lw);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
legend(h,{ 'untreated','peak signaling', 'adapted'},'Location','SouthEast')
pbaspect([3 2.5 1])

saveas(gcf, fullfile(dataDir, 'FRAPcombinedForSlide.png'));
saveas(gcf, fullfile(dataDir, 'FRAPcombinedForSlide.fig'));


%% single curve for slide

%untrIdx = [6 2; 11 4; 18 4; 19 4; 21 10];  
%peakIdx = [1 1; 3 2; 6 1; 7 1; 19 2; 21 1]; 
%adaptIdx = [3 4; 9 2; 19 3; 21 4; 22 3; 22 4]; % 19 3  

colors = lines(3);
idx = {[2 4], [1 1], [3 4]};
shapeIndices = [1 1 2];

shift = [5 0 -2];

% NEXT: change fit here to get nice fit for figure
clf
hold on
h = [];
for i = 2

    ei = idx{i}(1);
    fi = idx{i}(2);
    shapeIdx = shapeIndices(i);
    disp([FRAPdirs{ei} '/' results{ei}{fi}.description])

    tracesnorm = results{ei}{fi}.tracesNucNorm;%(shapeIdx,:);
    t = results{ei}{fi}.tres*(0:size(tracesnorm,2)-1);
    tmax = results{ei}{fi}.tmax;
    tlim = round(max(tmax)*results{ei}{fi}.tres);
    frapframe = results{ei}{fi}.frapframe;

    t = t + shift(i)*results{ei}{fi}.tres;
    lw = 1;
    h(i) = plot(t, tracesnorm(shapeIdx,:),...
                            'Color',colors(i,:),'LineWidth',lw)

    fitval = fitFRAP2(results{ei}{fi}); % change to 3 for recent version
    A = fitval.A(shapeIdx,1);
    k = fitval.k(shapeIdx,1);
    B = fitval.B(shapeIdx,1);

    func = @(p1,p2,p3,x) p1*(1-exp(-x*p2)) + p3;
    fitcurve = func(A,k,B,t(frapframe:end)-t(frapframe));  
    
    lw = 3;
    plot(t(frapframe:end),fitcurve,...
                            'Color', colors(i,:),'LineWidth',lw)
end

hold off

fs = 26;
yticks([0.1 0.3 0.5]);
%ylabel('relative recovery', 'FontSize',fs, 'FontWeight','Bold');
%xlabel('time (sec)', 'FontSize',fs, 'FontWeight','Bold');
xlim([0 1000]);
ylim([0 0.6]);
set(gcf,'color','w');
set(gca, 'LineWidth', lw);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
%legend(h,{ 'untreated','peak signaling', 'adapted'},'Location','SouthEast')
pbaspect([3 2.5 1])

saveas(gcf, fullfile(dataDir, 'FRAPsingleForSlide.png'));
saveas(gcf, fullfile(dataDir, 'FRAPsingleForSlide.fig'));

%% MODEL FIT actual values

% measured :
% recovery amplitudes A
% nuc:cyt ratios R
% inverse recovery time k
input = {};
output = {};

Aname = 'Acorr'; 
Acname = 'AcorrCyt';
% Aname = 'Araw';
% Acname = 'ArawCyt';
instruct = {untreated_in, peak_in, adapted_in};
instructp = {untreatedLMB_in, peakLMB_in, adaptedLMB_in};

disp('-------------------------------------------------');

for i = 1:3

    input{i} = struct();

    input{i}.A = mean(instruct{i}.(Aname));
    input{i}.Ap = mean(instructp{i}.(Aname));
    input{i}.R = mean(instruct{i}.ncrinit);
    input{i}.Rnb = mean(instruct{i}.ncrfin);
    input{i}.Rp = mean(instructp{i}.ncrinit);
    input{i}.k = mean(instruct{i}.k);
    input{i}.kp = mean(instructp{i}.k);
    
    input{i}.Ac = mean(instruct{i}.(Acname));
    input{i}.Acp = mean(instructp{i}.(Acname));
    
    % measured standard errors
    input{i}.sigA = std(instruct{i}.(Aname))/sqrt(sum(instruct{i}.Ncells)-1);
    input{i}.sigAp = std(instructp{i}.(Aname))/sqrt(sum(instructp{i}.Ncells)-1);
    input{i}.sigR = std(instruct{i}.ncrinit)/sqrt(sum(instruct{i}.Ncells)-1);
    input{i}.sigRp = std(instructp{i}.ncrinit)/sqrt(sum(instructp{i}.Ncells)-1);
    input{i}.sigk = std(instruct{i}.k)/sqrt(sum(instruct{i}.Ncells)-1);
    input{i}.sigkp = std(instructp{i}.k)/sqrt(sum(instructp{i}.Ncells)-1);
    
    input{i}.sigAc = std(instruct{i}.(Acname))/sqrt(sum(instruct{i}.Ncells)-1);
    input{i}.sigAcp = std(instructp{i}.(Acname))/sqrt(sum(instructp{i}.Ncells)-1);
end

alpha = 1.13;
measuredInput = input;
% 
% %%
% output = fitKineticModel4(input);%,alpha);
% output2 = fitKineticModel3(input,alpha);
% output2 = fitKineticModelFixKin(input);
% output2 = fitKineticModelFixKinNoCytSeq(input);

%% try fitting cs = 0 model to data without LMB

input = measuredInput;
%input{1}.Ac = input{1}.Ac*1.3;
%input{2}.A = 0.4;
%input{2}.A = 0.4;
%input{3}.A = 0.65;
input{3}.Ac = 0.17;
% ARE THESE ADJUSTMENTS PLAUSIBLE?
% CAN I JUST FIT TO ADJUSTED AMPLITUDES? THINK SO

% predicted Rrat should be adjusted for bn bc bleaching

output = fitKineticModelNoLMB(input);
% also try ns = 0, cs = free

for i = 1:3
    N = 2;
    disp('--------------------------');
    fprintf('var \tfit \t measure\n');
    disp('--------------------------');
    fprintf(['An:\t' num2str([An(output{i}.kin, output{i}.kout, output{i}.cs, output{i}.ns) input{i}.A],N) '\n']);
    fprintf(['Ac:\t' num2str([Ac(output{i}.kin, output{i}.kout, output{i}.cs, output{i}.ns) input{i}.Ac],N) '\n']);
    fprintf(['k:\t' num2str([k(output{i}.kin, output{i}.kout) input{i}.k],N) '\n']);
    
    % predicted R before bleach
    Rfit = R(output{i}.kin, output{i}.kout, output{i}.cs, output{i}.ns);
    % predicted R after bleach
    Rnbfit = Rnb(output{i}.kin, output{i}.kout, output{i}.cs, output{i}.ns);
    Rrat = Rnbfit/Rfit;
    fprintf(['R:\t' num2str([Rfit input{i}.R Rfit/input{i}.R],N)  '\n']);
    fprintf(['Rnb:\t' num2str([Rnbfit input{i}.Rnb Rnbfit/input{i}.Rnb],N) '\n']);
    fprintf(['Rrat:\t' num2str([Rrat input{i}.Rnb/input{i}.R],N) '\n']);
end

%% amplitude based model

manualInput = {};
manualInput{3} = struct('A',0.65, 'Ac',0.15,'k',0.004,'Ap',0.3,'Acp',0.3,'kp',0.001,...
                        'sigA',1, 'sigAc',1,'sigk',1,'sigAp',1,'sigAcp',1,'sigkp',1);
                    
manualInput{1} = manualInput{3};
manualInput{2} = manualInput{3};

input = measuredInput;
% input{3}.Ac = 0.17;
% input{3}.Ap = 0.1;
% input{3}.kp = 0.0005;
%input = manualInput;
output = fitKineticModelNew(input);

%output = fitKineticModelNew(manualInput);

% consistency check
kap = @(kin, kout) kin/(kin+kout);
An = @(kin, kout, cs, ns) kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs)/(ns + kap(kin,kout)*(1-ns-cs));
Ac = @(kin, kout, cs, ns) kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs)/(cs + (1-kap(kin,kout))*(1-ns-cs));
k = @(kin, kout) kin + kout;
R = @(kin, kout, cs, ns) (ns + kap(kin,kout)*(1-ns-cs))/(cs + (1-kap(kin,kout))*(1-ns-cs));
Rnb = @(kin, kout, cs, ns) (kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs))/(cs + (1-kap(kin,kout))^2*(1-ns-cs));

%Anb = @(kin, kout, cs, ns, bn, bc) kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs)/(ns + kap(kin,kout)*(1-ns-cs));
%Acb = @(kin, kout, cs, ns, bn, bc) kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs)/(cs + (1-kap(kin,kout))*(1-ns-cs));

for i = 1:3
    N = 2;
    disp('--------------------------');
    fprintf('var \tfit \t measure\n');
    disp('--------------------------');
    fprintf(['An:\t' num2str([An(output{i}.kin, output{i}.kout, output{i}.cs, output{i}.ns) input{i}.A],N) '\n']);
    fprintf(['Ac:\t' num2str([Ac(output{i}.kin, output{i}.kout, output{i}.cs, output{i}.ns) input{i}.Ac],N) '\n']);
    fprintf(['k:\t' num2str([k(output{i}.kin, output{i}.kout) input{i}.k],N) '\n']);
    fprintf(['Ap:\t' num2str([An(output{i}.kin, output{i}.koutp, output{i}.cs, output{i}.ns) input{i}.Ap],N) '\n']);
    fprintf(['Acp:\t' num2str([Ac(output{i}.kin, output{i}.koutp, output{i}.cs, output{i}.ns) input{i}.Acp],N) '\n']);
    fprintf(['kp:\t' num2str([k(output{i}.kin, output{i}.koutp) input{i}.kp],N) '\n']);
    fprintf(['R:\t' num2str([Rfit input{i}.R Rfit/input{i}.R],N)  '\n']);
    fprintf(['Rnb:\t' num2str([Rnbfit input{i}.Rnb Rnbfit/input{i}.Rnb],N) '\n']);
    fprintf(['Rrat:\t' num2str([Rrat input{i}.Rnb/input{i}.R],N) '\n']);
    fprintf(['Rp:\t' num2str([R(output{i}.kin, output{i}.koutp, output{i}.cs, output{i}.ns) input{i}.Rp],N) '\n']);
end

%% consistency check R-based model

% changing R
output = output2;
alpha = output2{1}.alpha;

kap = @(kin, kout) kin/(kin+kout);
A = @(kin, kout, cs, ns) kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs)/(ns + kap(kin,kout)*(1-ns-cs));
R = @(kin, kout, cs, ns) alpha*(kout*ns + kin*(1-cs))/(kin*cs + kout*(1-ns));
k = @(kin, kout) kin + kout;

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

