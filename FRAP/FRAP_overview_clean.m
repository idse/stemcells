clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;

FRAPmetadata_clean

%% combine nuclear bleach measurements for same condition  

disp('------------');
pf = '171217c';
load(fullfile(dataDir, ['FRAPresults' pf '.mat']));

untrIdx = [6 2; 11 4; 18 4; 19 4; 21 10]; 
peakIdx = [1 1; 3 2; 6 1; 7 1; 19 2; 21 1]; 
adaptIdx = [3 4; 9 2; 19 3; 21 4; 22 3; 22 4]; 

untrIdxLMB = [2 3; 9 6; 9 7; 11 2; 11 3; 16 3];
peakIdxLMB = [9 1; 10 2; 11 1; 14 1; 21 5;];
adaptIdxLMB = [9 3; 10 1; 17 1; 21 8; 23 5; 23 6]; 

allIdx = {untrIdx, peakIdx, adaptIdx, untrIdxLMB, peakIdxLMB, adaptIdxLMB};
allOut = {};
allLabels = {'untreated', 'peak', 'adapted',...
                    'untreatedLMB', 'peakLMB', 'adaptedLMB'};

debugplots = false;

% definitions of amplitudes in terms of intensities
ANm = @(Ni, Nb, Nf, B) (Nf - Nb)./(Ni - B);
BNm = @(Ni, Nb, Nf, B) (Nb - B)./(Ni - B);
ACm = @(Ci, Cb, Cf, B) (Cb - Cf)./(Ci - B);
BCm = @(Ci, Cb, Cf, B) (Cf - B)./(Ci - B);
bc  = @(Ci, Cb, Cf, B) (Cb - B)./(Ci - B);

% optional transformations accounting for microscope
% not used in the end for the sake of simplicity
Bcam = 110;
Bm = 40;
a = 1; b = (1-a); % a = 0.45 = optimal
fn = @(N, B) N - B;
fc = @(C, B) (C - B - b*Bm)/a;

for i = 1:numel(allIdx)
    
    legendstr = FRAPdirs(allIdx{i}(:,1));
    allOut{i} = makeInStruct(results, allIdx{i}, debugplots, legendstr, dataDir, allLabels{i},bleachCorrect);
    
    % definitions of measured intensities
    B = allOut{i}.bg;
    %Bcam = B; % switch between fixed and variable background subtraction
    Ni = fn(allOut{i}.Ninit,    Bcam);
    Nb = fn(allOut{i}.Nbleach,  Bcam);
    Nf = fn(allOut{i}.Nfin,     Bcam);
    Ci = fc(allOut{i}.Cinit,    Bcam);
    Cb = fc(allOut{i}.Cbleach,  Bcam);
    Cf = fc(allOut{i}.Cfin,     Bcam);

    bg = allOut{i}.bg;
    
    allOut{i}.Araw = (Nf-Nb)./Ni;
    allOut{i}.Braw = Nb./Ni;
    allOut{i}.bn = allOut{i}.Braw;
    
    allOut{i}.ArawCyt = (Cb-Cf)./Ci;
    allOut{i}.BrawCyt = Cf./Ci;
    allOut{i}.bc = Cb./Ci;

    allOut{i}.Acorr = allOut{i}.Araw./(allOut{i}.bc - allOut{i}.bn);
    allOut{i}.AcorrCyt = allOut{i}.ArawCyt./(allOut{i}.bc - allOut{i}.bn);

    allOut{i}.Rinit = Ni./Ci;
    allOut{i}.Rfin = Nf./Cf;
    allOut{i}.Rrat = allOut{i}.Rfin./allOut{i}.Rinit;
    
    % raw is tilde, corr is no tilde
    allOut{i}.RratCorr = allOut{i}.Acorr./(1 - allOut{i}.AcorrCyt); 
    
    allOut{i}.tau = 1./(60*allOut{i}.k);
end

untreated_in = allOut{1};
peak_in = allOut{2}; 
adapted_in = allOut{3}; 

untreatedLMB_in = allOut{4};
peakLMB_in = allOut{5};
adaptedLMB_in = allOut{6};

% for i = 1:3
%     disp(['-' num2str(i)])
%     disp(mean([allOut{i}.Rrat allOut{i}.RratCorr]))
%     disp(std([allOut{i}.Rrat allOut{i}.RratCorr]))
% end

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

for i = 1:3

    input{i} = struct();

    input{i}.A = mean(instruct{i}.(Aname));
    input{i}.r = mean(instruct{i}.Rinit);
    input{i}.rp = mean(instruct{i}.Rfin); %uncorrected
    %input{i}.rp = mean(instruct{i}.ncrinit.*instruct{i}.(Aname)./(1-instruct{i}.(Acname))); %corrected
    
    input{i}.k = mean(instruct{i}.k);
    input{i}.N = sum(instruct{i}.Ncells);
    
    input{i}.Ac = mean(instruct{i}.(Acname));
    input{i}.Aclmb = mean(instructp{i}.(Acname));

    % measured standard errors
    input{i}.sigA = std(instruct{i}.(Aname))/sqrt(sum(instruct{i}.Ncells)-1);
    input{i}.sigr = std(instruct{i}.Rinit)/sqrt(sum(instruct{i}.Ncells)-1);
    input{i}.sigk = std(instruct{i}.k)/sqrt(sum(instruct{i}.Ncells)-1);
    
    input{i}.sigAc = std(instruct{i}.(Acname))/sqrt(sum(instruct{i}.Ncells)-1);
end

measuredInput = input;

kin = [1 1 1]*0.00089;
[output,res] = fitKineticModelFixKinFitAlpha(input, kin);
di = 4; % equations per condition

ressq = res.^2;

% consistency check
kap = @(kin, kout) kin/(kin+kout);
An = @(kin, kout, cs, ns) kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs)/(ns + kap(kin,kout)*(1-ns-cs));
Ac = @(kin, kout, cs, ns) kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs)/(cs + (1-kap(kin,kout))*(1-ns-cs));
k = @(kin, kout) kin + kout;
R = @(kin, kout, cs, ns) (ns + kap(kin,kout)*(1-ns-cs))/(cs + (1-kap(kin,kout))*(1-ns-cs));
Rnb = @(kin, kout, cs, ns) (kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs))/(cs + (1-kap(kin,kout))^2*(1-ns-cs));
Rcb = @(kin, kout, cs, ns) (ns + kap(kin,kout)^2*(1-ns-cs))/((1-kap(kin,kout))*kap(kin,kout)*(1-ns-cs));

%Anb = @(kin, kout, cs, ns, bn, bc) kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs)/(ns + kap(kin,kout)*(1-ns-cs));
%Acb = @(kin, kout, cs, ns, bn, bc) kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs)/(cs + (1-kap(kin,kout))*(1-ns-cs));

disp('----------------------------------------------------');
fprintf('var \tfit \t measure \terror\tdiff\tressq\n');

for i = 1:3
    if ~isfield(output{i},'beta')
        output{i}.beta = 1;
    end
    N = 2;
    disp('----------------------------------------------------');
    AnFit = An(output{i}.kin, output{i}.kout, output{i}.cs, output{i}.ns);
    AcFit = output{i}.beta*Ac(output{i}.kin, output{i}.kout, output{i}.cs, output{i}.ns);
    kfit = k(output{i}.kin, output{i}.kout);
    fprintf(['An:\t' num2str([AnFit input{i}.A],N) '\t\t' num2str([input{i}.sigA  AnFit-input{i}.A ressq(1 + di*(i-1))],N) '\n']);
    fprintf(['Ac:\t' num2str([AcFit input{i}.Ac],N) '\t\t' num2str([input{i}.sigAc AcFit-input{i}.Ac ressq(2 + di*(i-1))],N) '\n']);
    fprintf(['k:\t' num2str([kfit input{i}.k],N) '\t\t' num2str([input{i}.sigk kfit-input{i}.k ressq(3 + di*(i-1))],N) '\n']);
    fprintf(['kappa:\t' num2str([kap(output{i}.kin, output{i}.kout)],N) '\n']);
    
    % predicted R before bleach
    Rfit = R(output{i}.kin, output{i}.kout, output{i}.cs, output{i}.ns);
    % predicted R after nuclear bleach
    Rpfit = Rnb(output{i}.kin, output{i}.kout, output{i}.cs, output{i}.ns);
    % predicted R after cyto bleach
    Rcbfit = Rcb(output{i}.kin, output{i}.kout, output{i}.cs, output{i}.ns);

    Rrat = Rpfit/Rfit;
    fprintf(['R:\t' num2str([Rfit input{i}.r/output{i}.alpha],N)  '\t\t' num2str([input{i}.sigr Rfit-input{i}.r/output{i}.alpha ressq(4 + di*(i-1))],N) '\n']);

    disp('-consistency-');
    fprintf(['Ac/An:\t' num2str([AcFit/AnFit input{i}.Ac/input{i}.A],N)  '\n']);
    %fprintf(['rp:\t' num2str([Rpfit*output{i}.alpha input{i}.rp],N) '\n']);

    %fprintf(['Rrat:\t' num2str([Rrat input{i}.rp/input{i}.r],N) '\n']);
    % corrected version:
    fprintf(['Rrat:\t' num2str([Rrat input{i}.A./(1-input{i}.Ac/output{i}.beta)],N) '\n']);
end

%% sensitivity matrices for paper

num2str(output{1}.sensitivity(1:4,1:3),2)
num2str(output{1}.sensitivity(5:8,4:6),2)
num2str(output{1}.sensitivity(9:12,7:9),2)

%% p - values for paper
% e.g. https://www.statsdirect.co.uk/help/parametric_methods/utt.htm

i = 2;
j = 3;

varname = 'kout';

mu1 = output{i}.(varname);
mu2 = output{j}.(varname);

n1 = input{i}.N;
n2 = input{j}.N;

s1 = sqrt(n1-1)*output{i}.(['sig' varname]);
s2 = sqrt(n1-1)*output{j}.(['sig' varname]);

a = s1^2/n1;
b = s2^2/n2;
t = abs(mu1 - mu2)/sqrt(a+b);
nu = (a+b)^2/(a^2/(n1-1) + b^2/(n2-1));
p = 1 - tcdf(t,nu)

% %%
% 
% [h,p] = ttest2(allOut{i}.Acorr, allOut{j}.Acorr);
% 
% mu1 = mean(allOut{j}.Acorr);
% mu2 = mean(allOut{j}.Acorr);
% s1 = std(allOut{i}.Acorr);
% s2 = std(allOut{j}.Acorr);
% n1 = allOut{i}.N;
% n2 = allOut{j}.N;
% t = sqrt(n2)*(mu2-mu1)/s2;
% 1-tcdf(t, n2-1)

%% box plots for R'/R just nuclear

figure, 
hold on
valsvec = [];
idxvec = [];
coloridx = [];
color = lines(2);

lw = 3;
fs = 28;
whiskerlength = Inf;
for i = 1:3
    valsvec = cat(1,valsvec, allOut{i}.RratCorr);
    coloridx = cat(1,coloridx, allOut{i}.RratCorr*0+1);
    idxvec = cat(1,idxvec, allOut{i}.RratCorr*0+2*i-1);
end 
h = boxplot(valsvec, idxvec,'notch','off','Widths',0.5,...
                        'ColorGroup',coloridx,'Colors',color(2,:),...
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
ylabel("R'/R");
box off
pbaspect([3 2 1]);
axis square
hold off
%saveas(gcf, fullfile(dataDir,['rcompareJustNuclear' titlestr '_boxplot.png']));

%% plot parameter values

whiskerlength = Inf;

disp('display measured parameters');

measured = {'A','Acorr','k','x'};% 'rinit','rfin', % w and w/o LMB
legloc = {'NorthWest','NorthEast','NorthEast','NorthWest'};
yranges = {[0 0.8],[0 1],[0 0.0085]};%[0 20],[0 5]};
titles = {'recovery amplitude A','recovery amplitude A','recovery rate k'};
for j = 2:3

    name = measured{j};
    
    vals = {untreated_in.(name), peak_in.(name), adapted_in.(name)}; 
    means = [mean(untreated_in.(name)) mean(peak_in.(name)) mean(adapted_in.(name))]; 
    errorvals = 2*[   std(untreated_in.(name))/sqrt(sum(untreated_in.Ncells)-1)...
                    std(peak_in.(name))/sqrt(sum(peak_in.Ncells)-1)...
                    std(adapted_in.(name))/sqrt(sum(adapted_in.Ncells)-1)]; 
                %[std(untreated_in.(name)) std(peak_in.(name)) std(adapted_in.(name))]; 
    
    % box plots
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
    % following line for single color 
    coloridx = 0*coloridx+1;
    h = boxplot(valsvec, idxvec,'notch','off','Widths',0.5,...
                            'ColorGroup',coloridx,'Colors',color(2,:),...
                            'Whisker',whiskerlength); % 1.5 default
    set(h,'linew',lw);
    set(gca, 'LineWidth', lw);
    set(gcf,'color','w');
    if ~isempty(yranges{j})
        ylim(yranges{j});
    end
    ylabel(titles{j});
    xticks([1 2 3]);
    set(gca,'XTickLabel', {'untreat', 'peak', 'adapt'});
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    box off
    axis square
    %saveas(gcf, fullfile(dataDir,[name '_boxplotreduced.png']));
end

%% display amounts of input data

Nmovies = [untreated_in.Nmovies peak_in.Nmovies adapted_in.Nmovies]
NmoviesLMB = [untreatedLMB_in.Nmovies peakLMB_in.Nmovies adaptedLMB_in.Nmovies]

Ncells = [sum(untreated_in.Ncells) sum(peak_in.Ncells) sum(adapted_in.Ncells)]
NcellsLMB = [sum(untreatedLMB_in.Ncells) sum(peakLMB_in.Ncells) sum(adaptedLMB_in.Ncells)]

%% overlay some FRAP curves

colors = lines(3);
idx = {[2 4], [1 1], [3 4]};
shapeIndices = [1 1 2];

shift = [5 0 -2]; %shift t=0 to the same place

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

%saveas(gcf, fullfile(dataDir, 'FRAPcombinedForSlide.png'));
%saveas(gcf, fullfile(dataDir, 'FRAPcombinedForSlide.fig'));


%% single curve for figure

colors = lines(3);
idx = {[2 4], [1 1], [3 4]};
shapeIndices = [1 1 2];

shift = [5 0 -2];

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
                            'Color',colors(2,:),'LineWidth',lw)

    fitval = fitFRAP2(results{ei}{fi}); % change to 3 for recent version
    A = fitval.A(shapeIdx,1);
    k = fitval.k(shapeIdx,1);
    B = fitval.B(shapeIdx,1);

    func = @(p1,p2,p3,x) p1*(1-exp(-x*p2)) + p3;
    fitcurve = func(A,k,B,t(frapframe:end)-t(frapframe));  
    
    lw = 3;
    plot(t(frapframe:end),fitcurve,...
                            'Color', colors(2,:),'LineWidth',lw)
end

hold off

fs = 28;
yticks([0.1 0.3 0.5]);
ylabel('relative recovery', 'FontSize',fs, 'FontWeight','Bold');
xlabel('time (sec)', 'FontSize',fs, 'FontWeight','Bold');
xlim([0 1200]);
ylim([0 0.5]);
set(gcf,'color','w');
set(gca, 'LineWidth', lw);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
%legend(h,{ 'untreated','peak signaling', 'adapted'},'Location','SouthEast')
%pbaspect([3 2.5 1])
axis square

%saveas(gcf, fullfile(dataDir, 'FRAPsingleForSlide.png'));
%saveas(gcf, fullfile(dataDir, 'FRAPsingleForSlide.fig'));

%% visualize model fit 

prefix = ['kinfixed' num2str(num2str(kin(1))) '_'];
outputcomb = cat(2,output{:});
inferred = {'kin','kout','cs','ns'};
vals = {};
errorvals = {};
titlestr = {'nuclear exchange','sequestration'};
legendpos = {'NorthEast','NorthWest'};
legendstr = {{'export','import'},{'nuclear','cytoplasm'}};

% % for relative changes, see below too
% ylabelstr = {'rate (10^{-3} sec^{-1})', '% change'};
% ylimval = {[0 7],[-100 100]};

% for absolute changes
ylabelstr = {'rate (10^{-3} sec^{-1})', 'Smad4 fraction'};
ylimval = {[0 7],[0 0.6]};

for j = 2
    
    figure,
    hold on
    %s = 10^3;
    s={10^3,1};
    k=1;
    for i = 2*j-1:2*j
        name = inferred{i};
        errname = ['sig' name];
        vals{k} = s{j}*[outputcomb.(name)];
        errorvals{k} = s{j}*[outputcomb.(errname)];
        k = k+1;
    end
%     % this is to get relative changes in sequestration
%     % for absolute, comment out
%     if j == 2
%         for ki = 1:k-1
%             errorvals{ki} = 100*(vals{ki}/vals{ki}(1)).*sqrt((errorvals{ki}(1)./vals{ki}(1))^2 + (errorvals{ki}./vals{ki}).^2);
%             errorvals{ki}(1) = 0;
% %             errorvals{ki} = 100*errorvals{ki}./vals{ki}(1);
%             vals{ki} = 100*(vals{ki}./vals{ki}(1) - 1);
%         end
%     end
    hold off

    fs = 26;
    colors = lines(2);
    hold on
    bins = (1:3);
    lw = 3;
    
    if j == 1
        % for the kin fixed model
        b = bar(vals{2},0.9);%'FaceColor',[0.1 0.5 0.1]);
        b(1).FaceColor = colors(1,:);
        errorbar(bins, vals{2}, errorvals{2},...
                            'k','linestyle','none','linewidth',lw);
    else
        b = bar(cat(1,vals{2}, vals{1})',0.9);%'FaceColor',[0.1 0.5 0.1]);
        b(1).FaceColor = colors(1,:);
        b(2).FaceColor = colors(2,:);
        errorbar(bins - 0.15, vals{2}, errorvals{2},...
                            'k','linestyle','none','linewidth',lw);
        errorbar(bins + 0.15, vals{1}, errorvals{1},...
                            'k','linestyle','none','linewidth',lw);
        legend(legendstr{j},'Location',legendpos{j});
    end

    hold off
    set(gcf,'color','w');
    set(gca, 'LineWidth', lw);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    box off
    set(gca,'XTickLabel', {'untreated', 'peak', 'adapt'});

    ylabel(ylabelstr{j});
    %title(titlestr{j})
    fname = [prefix titlestr{j} '_valuesRel2.png'];

    ylim(ylimval{j});
    xlim([0.5 3.5]);
    ylim(ylimval{j});

    axis square
    %saveas(gcf, fullfile(dataDir, fname));
end
