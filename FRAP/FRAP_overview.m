clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
dataDir = '/Users/idse/data_tmp/0_FRAP/';
dataDir = '/Volumes/IdseData4/0_FRAP/';

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

results{6}{1}.good = [1 0];
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

%save(fullfile(dataDir, ['FRAPresults' pf '.mat']), 'results');


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
    allOutCyt{i}.AcorrCyt = allOutCyt{i}.ArawCyt./(allOutCyt{i}.bn - allOutCyt{i}.bc);
    
    allOutCyt{i}.Rfin = (allOutCyt{i}.Nfin-bg)./(allOutCyt{i}.Cfin-bg);
    allOutCyt{i}.Rinit = (allOutCyt{i}.Ninit-bg)./(allOutCyt{i}.Cinit-bg);
    allOutCyt{i}.Rrat = allOutCyt{i}.Rfin./allOutCyt{i}.Rinit;
    
    allOutCyt{i}.RratCorr = (1-allOutCyt{i}.Acorr)./allOutCyt{i}.AcorrCyt;
end

%
% 
% %%
% i = 1
% r = (allOutCyt{i}.Cfin - allOutCyt{i}.Cbleach)./(allOutCyt{i}.Cinit-allOutCyt{i}.bg);
% [r allOutCyt{i}.Cfin allOutCyt{i}.Cbleach allOutCyt{i}.Cinit allOutCyt{i}.bg]

%% combine nuclear bleach measurements for same condition  

disp('------------');
pf = '171217c';
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
%peakIdxLMB = [9 1; 10 2; 11 2; 13 1; 14 1; 19 1; 21 5; 21 6];
untrIdxLMB = [2 3; 9 6; 9 7; 11 2; 11 3; 16 3];
peakIdxLMB = [9 1; 10 2; 11 1; 14 1; 21 5;];
% 13 1; garbage trace: too much movement?
% 19 1; LMB not long enough, recovery too fast, bad cyto
%  21 6
adaptIdxLMB = [9 3; 10 1; 17 1; 21 8; 23 5; 23 6]; %16 2; 

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

%Ni = fn(allOut{i}.Ninit,    Bcam);

for i = 1:numel(allIdx)
    
    legendstr = FRAPdirs(allIdx{i}(:,1));
    allOut{i} = makeInStruct(results, allIdx{i}, debugplots, legendstr, dataDir, allLabels{i},bleachCorrect);
    
    % definitions of measured intensities
    B = allOut{i}.bg;
    %Bcam = B;
    Ni = fn(allOut{i}.Ninit,    Bcam);
    Nb = fn(allOut{i}.Nbleach,  Bcam);
    Nf = fn(allOut{i}.Nfin,     Bcam);
    Ci = fc(allOut{i}.Cinit,    Bcam);
    Cb = fc(allOut{i}.Cbleach,  Bcam);
    Cf = fc(allOut{i}.Cfin,     Bcam);

    bg = allOut{i}.bg;
    
%     allOut{i}.Araw = (allOut{i}.Nfin - allOut{i}.Nbleach)./(allOut{i}.Ninit-bg);
%     allOut{i}.Braw = (allOut{i}.Nbleach - bg)./(allOut{i}.Ninit-bg);
    allOut{i}.Araw = (Nf-Nb)./Ni;
    allOut{i}.Braw = Nb./Ni;
    allOut{i}.bn = allOut{i}.Braw;
    
%     allOut{i}.ArawCyt = (allOut{i}.Cbleach-allOut{i}.Cfin)./(allOut{i}.Cinit-bg);
%     allOut{i}.BrawCyt = (allOut{i}.Cfin - bg)./(allOut{i}.Cinit-bg);
%     allOut{i}.bc = (allOut{i}.Cbleach - bg)./(allOut{i}.Cinit-bg);
    allOut{i}.ArawCyt = (Cb-Cf)./Ci;
    allOut{i}.BrawCyt = Cf./Ci;
    allOut{i}.bc = Cb./Ci;

    allOut{i}.Acorr = allOut{i}.Araw./(allOut{i}.bc - allOut{i}.bn);
    allOut{i}.AcorrCyt = allOut{i}.ArawCyt./(allOut{i}.bc - allOut{i}.bn);

    %allOut{i}.Rfin = (allOut{i}.Nfin-bg)./(allOut{i}.Cfin-bg);
    %allOut{i}.Rinit = (allOut{i}.Ninit-bg)./(allOut{i}.Cinit-bg);
    %allOut{i}.Rrat = allOut{i}.ncrfin./allOut{i}.ncrinit;
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

disp('-------------------------------------------------');
disp(['a = ' num2str(a) ', b = ' num2str(b)]);

for i = 1:3

    input{i} = struct();

    input{i}.A = mean(instruct{i}.(Aname));
    input{i}.Almb = mean(instructp{i}.(Aname));
    input{i}.r = mean(instruct{i}.Rinit);
    
    input{i}.rp = mean(instruct{i}.Rfin); %uncorrected
    %input{i}.rp = mean(instruct{i}.ncrinit.*instruct{i}.(Aname)./(1-instruct{i}.(Acname))); %corrected
    
    input{i}.rlmb = mean(instructp{i}.Rinit);
    input{i}.k = mean(instruct{i}.k);
    input{i}.klmb = mean(instructp{i}.k);
    input{i}.N = sum(instruct{i}.Ncells);
    
    input{i}.Ac = mean(instruct{i}.(Acname));
    input{i}.Aclmb = mean(instructp{i}.(Acname));

    % measured standard errors
    input{i}.sigA = std(instruct{i}.(Aname))/sqrt(sum(instruct{i}.Ncells)-1);
    input{i}.sigAlmb = std(instructp{i}.(Aname))/sqrt(sum(instructp{i}.Ncells)-1);
    input{i}.sigr = std(instruct{i}.Rinit)/sqrt(sum(instruct{i}.Ncells)-1);
    input{i}.sigrlmb = std(instructp{i}.Rinit)/sqrt(sum(instructp{i}.Ncells)-1);
    input{i}.sigk = std(instruct{i}.k)/sqrt(sum(instruct{i}.Ncells)-1);
    input{i}.sigklmb = std(instructp{i}.k)/sqrt(sum(instructp{i}.Ncells)-1);
    
    input{i}.sigAc = std(instruct{i}.(Acname))/sqrt(sum(instruct{i}.Ncells)-1);
    input{i}.sigAclmb = std(instructp{i}.(Acname))/sqrt(sum(instructp{i}.Ncells)-1);
end

measuredInput = input;

% fit full model using LMB -- is r read out with same cutoff?

% manualInput = {};
% manualInput{3} = struct('A',0.65, 'Ac',0.15,'k',0.004,'Ap',0.3,'Acp',0.3,'kp',0.001,...
%                         'sigA',1, 'sigAc',1,'sigk',1,'sigAp',1,'sigAcp',1,'sigkp',1);
%                     
% manualInput{1} = manualInput{3};
% manualInput{2} = manualInput{3};

input = measuredInput;
%input{1}.Ac = 0.1;

% adjustments to play with
s = 1;
sklmb = 1;
sAlmb = 1;
sAclmb = 1;
for i = 1:3
    input{i}.Ac = s*input{i}.Ac;
    input{i}.klmb = sklmb*input{i}.klmb;
    input{i}.Almb = sAlmb*input{i}.Almb;
    input{i}.Aclmb = sAclmb*input{i}.Aclmb;
end
%input{1}.A = 0.5;
%input{2}.Ac = s*input{2}.Ac;
%input{2}.Almb= 0.05;
%input{2}.Almb =0.1;
%input{2}.Aclmb =0.6;
%input{2}.Aclmb = s*input{2}.Aclmb;
%input{2}.klmb = 0.001;
%input{2}.rlmb = 8;
% 
% input{3}.Almb= 0.3;
% input{3}.Aclmb = 0.7;
% input{3}.klmb = 0.001;
% input{3}.A =0.8;
%input{3}.Ac = 0.3;

% FIT LMB DATA WITH REDUCED MODEL
% FIT LMB DATA WITH FIXED KOUTP

%[output,res] = fitKineticModelNewFitAlphaFixKoutp(input);
% [output,res] = fitKineticModelNewFitAlphaFixKin(input);

% 180516 Note: see if assuming different cs for LMB & no LMB can give a
% good fit to support the conclusion that LMB disturbs equilibrium of cs?
%[output,res] = fitKineticModelNewFitAlphaCsPrime(input);

[output,res] = fitKineticModelNewFitAlpha(input);
reduced = false;
di = 8; % equations per condition

% [output, res] = fitKineticModelFitAlphaReduced(input);
% di = 6;

%output = fitKineticModelFitAlphaBeta(input);
% output = fitKineticModelNew(input);

% % reduced models
cs = 0.4*[1 1 1];
ns = 0.03*[1 1 1];
kin = [1 1 1]*0.00089;
%[output,res] = fitKineticModelFixKinFitAlpha(input, kin);
%[output,res] = fitKineticModelFixCsFitAlpha(input, cs);
% [output,res] = fitKineticModelFixNsFitAlpha(input, ns);
% reduced = true;
% di = 4; % equations per condition

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
    %fprintf(['kaplmb:\t' num2str([kap(output{i}.kin, output{i}.koutp)],N) '\n']);
    
    if ~reduced
        AnlmbFit = An(output{i}.kin, output{i}.koutp, output{i}.cs, output{i}.ns);
        AclmbFit = output{i}.beta*Ac(output{i}.kin, output{i}.koutp, output{i}.cs, output{i}.ns);
        fprintf(['Almb:\t' num2str([AnlmbFit input{i}.Almb],N) '\n']);
        fprintf(['Aclmb:\t' num2str([AclmbFit input{i}.Aclmb],N) '\n']);
        fprintf(['klmb:\t' num2str([k(output{i}.kin, output{i}.koutp) input{i}.klmb],N) '\n']);
    end
    
    % predicted R before bleach
    Rfit = R(output{i}.kin, output{i}.kout, output{i}.cs, output{i}.ns);
    % predicted R after nuclear bleach
    Rpfit = Rnb(output{i}.kin, output{i}.kout, output{i}.cs, output{i}.ns);
    % predicted R after cyto bleach
    Rcbfit = Rcb(output{i}.kin, output{i}.kout, output{i}.cs, output{i}.ns);

    Rrat = Rpfit/Rfit;
    fprintf(['R:\t' num2str([Rfit input{i}.r/output{i}.alpha],N)  '\t\t' num2str([input{i}.sigr Rfit-input{i}.r/output{i}.alpha ressq(4 + di*(i-1))],N) '\n']);
    if ~reduced
        RfitLMB = R(output{i}.kin, output{i}.koutp, output{i}.cs, output{i}.ns);
        fprintf(['Rlmb:\t' num2str([RfitLMB input{i}.rlmb/output{i}.alpha],N)  '\t\t' num2str([input{i}.sigrlmb RfitLMB-input{i}.rlmb/output{i}.alpha ressq(di*i)],N) '\n']);
    end
    
    disp('-consistency-');
    fprintf(['Ac/An:\t' num2str([AcFit/AnFit input{i}.Ac/input{i}.A],N)  '\n']);
    if ~reduced
        fprintf(['Ac/AnL:\t' num2str([AclmbFit/AnlmbFit input{i}.Aclmb/input{i}.Almb],N)  '\n']);
    end
    %fprintf(['rp:\t' num2str([Rpfit*output{i}.alpha input{i}.rp],N) '\n']);

    % no point in displaying the corrected version of R', 
    % it is just multiplying r and Rrat
    %fprintf(['rpcorr:\t' num2str([Rpfit*output{i}.alpha input{i}.r*input{i}.A./(1-input{i}.Ac)],N) '\n']);
    
    %fprintf(['Rrat:\t' num2str([Rrat input{i}.rp/input{i}.r],N) '\n']);
    % corrected version:
    fprintf(['Rrat:\t' num2str([Rrat input{i}.A./(1-input{i}.Ac/output{i}.beta)],N) '\n']);
    %fprintf(['Rcb:\t' num2str([Rcbfit*output{i}.alpha mean(allOutCyt{i}.ncrinit)],N)...
    %    '(' num2str(std(allOutCyt{i}.ncrinit),N) ')\n']);
    Nc = numel(allOutCyt{i}.RratCorr);
    fprintf(['Rcbrat:\t' num2str([Rcbfit/Rfit...
        mean(allOutCyt{i}.RratCorr)],N) '(' num2str(std(allOutCyt{i}.RratCorr)/sqrt(Nc),N) ')\n']);
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

%%
disp('-----------------------------');
fprintf(['Cinit\tCb \tCfin \tB\t Ac\n']);
for i = 1:6
    AcTest = (allOut{i}.Cbleach-allOut{i}.Cfin)./(allOut{i}.Cinit - allOut{i}.bg);
    AcTestMean = (mean(allOut{i}.Cbleach)-mean(allOut{i}.Cfin))./(mean(allOut{i}.Cinit) - mean(allOut{i}.bg));
    disp(num2str([mean([allOut{i}.Cinit allOut{i}.Cbleach allOut{i}.Cfin allOut{i}.bg AcTest])],'%.2f\t'))
end


%% display relevant values for cleanup

names = {'Araw','Braw','ArawCyt','BrawCyt','bc','Acorr','AcorrCyt','k',...
            'Rinit','ncrfin','Rrat','tau'};
fs = 15;
for i = 8
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

%% display relevant values for cleanup cytoplasm

% k is nuclear recovery time??
% A is the nuclear amplitude

names = {'Araw','Braw','ArawCyt','BrawCyt','bc','Acorr','AcorrCyt','k',...
            'Rinit','Rfin','Rrat','tau'};
fs = 15;
for i = 8%[3 7]%[1 3 8]%9:10%1:8
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

%%
disp('compare nuclear:cytoplasmic recovery');

titlestr = 'rcompare';

figure,
hold on
    
vals = {};
vals2 = {};

disp("Rinit_CBmean Rinit_NBmean Rfin_CBmean Rfin_NBmean Rrat_CBmean Rrat_NBmean Rrat_NBmeanCorr]");

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
    
    vals{i} = allOut{i}.RratCorr;
    vals2{i} = allOutCyt{i}.RratCorr;
    
    %[allOut{i}.ncrinit allOut{i}.ncrfin allOut{i}.ncrfin./allOut{i}.ncrinit]
    Rinit_NBmean = nanmean(allOut{i}.ncrinit);
    Rfin_NBmean = nanmean(allOut{i}.ncrfin);
    Rrat_NBmean = nanmean(allOut{i}.Rrat);
    Rrat_NBmeanCorr = nanmean(allOut{i}.Acorr./(1-allOut{i}.AcorrCyt));
    
    scatter(allOut{i}.Rrat*0 + i, allOut{i}.Rrat);
    %scatter(allOutCyt{i}.Rrat*0 + i+0.5, allOutCyt{i}.RratCorr);
    disp(num2str([Rinit_CBmean Rinit_NBmean Rfin_CBmean Rfin_NBmean Rrat_CBmean Rrat_NBmean Rrat_NBmeanCorr],2))
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
%saveas(gcf, fullfile(dataDir,['rcompare' titlestr '_boxplot.png']));

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
    valsvec = cat(1,valsvec, vals{i});
    coloridx = cat(1,coloridx, vals{i}*0+1);
    idxvec = cat(1,idxvec, vals{i}*0+2*i-1);
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
saveas(gcf, fullfile(dataDir,['rcompareJustNuclear' titlestr '_boxplot.png']));



%% plot parameter values

whiskerlength = Inf;

disp('display measured parameters');

measured = {'A','Acorr','k','x'};% 'rinit','rfin', % w and w/o LMB
legloc = {'NorthWest','NorthEast','NorthEast','NorthWest'};
yranges = {[0 0.8],[0 1],[0 0.0085]};%[0 20],[0 5]};
titles = {'recovery amplitude A','recovery amplitude A','recovery rate k'};
for j = 3

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

    figure, 
    colors = lines(2);
    fs = 28;
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
    saveas(gcf, fullfile(dataDir,[name '_boxplotreduced.png']));
end


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
    k = 1;
    
    for j = 1:size(condIdx,1)
        
        idx = condIdx(j,:);
        oibfile = oibfiles{idx(1)}{idx(2)};
        ngood = sum(results{idx(1)}{idx(2)}.good);
        fprintf([num2str(idx(1)) ' ' num2str(idx(2)) '\t' oibfile '(' num2str(ngood) ') \t\t\t(' FRAPdirs{idx(1)} ')\n']);
        
        krange = k:k+ngood-1;
        disp(num2str([allOut{i}.Ninit(krange,:) allOut{i}.Nbleach(krange,:) allOut{i}.Nfin(krange,:)...
            allOut{i}.Cinit(krange,:) allOut{i}.Cbleach(krange,:) allOut{i}.Cfin(krange,:) allOut{i}.bg(krange,:)],4))
        k = k+ngood;
    end
    %disp('-');
    %fprintf('Ninit \t Nbleach \t Nfin \t Cinit \t Cbleach \t Cfin \t bg\n');
    
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
    legend({'nuclear (N_s)', 'cytoplasmic (C_s)'},'Location','NorthEast','FontSize',20)
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

saveas(gcf, fullfile(dataDir, 'FRAPsingleForSlide.png'));
saveas(gcf, fullfile(dataDir, 'FRAPsingleForSlide.fig'));

%% try fitting different models holding some parameter fixed
% to data without LMB

%% more consistency of fixed cs

csvals = 0:0.05:1;
nsvals = [];
kinvals = [];
koutvals = [];

for k = 1:numel(csvals)
    cs = csvals(k)*[1 1 1];
    output{k} = fitKineticModelFixCsFitAlpha(input, cs);
    for i = 1:3
        nsvals(k,i) = output{k}{i}.ns;
        kinvals(k,i) = output{k}{i}.kin;
        koutvals(k,i) = output{k}{i}.kout;
    end
end

good = all(koutvals > 0 & kinvals > 0 & (nsvals > 0 & nsvals < 1),2);
csvals(good)

%%
plot(csvals,nsvals.*good)
legend({'1','2','3'})
ylim([0 1])

%%
plot(csvals,kinvals.*good)
legend({'1','2','3'})
%ylim([0 1]);

%%
plot(csvals,koutvals)
legend({'1','2','3'})
%ylim([0 0.01]);

%% more consistency of fixed ns

nsvals = 0:0.001:0.05;
csvals = [];
kinvals = [];
koutvals = [];

for k = 1:numel(nsvals)
    ns = nsvals(k)*[1 1 1];
    output{k} = fitKineticModelFixNsFitAlpha(input, ns);
    for i = 1:3
        csvals(k,i) = output{k}{i}.cs;
        kinvals(k,i) = output{k}{i}.kin;
        koutvals(k,i) = output{k}{i}.kout;
    end
end

good = all(koutvals > 0 & kinvals > 0 & (csvals > 0 & csvals < 1),2);
nsvals(good)

%%
plot(nsvals,csvals.*good)
legend({'1','2','3'})
ylim([0 1])

%%
plot(nsvals,kinvals.*good)
legend({'1','2','3'})
%ylim([0 1]);

%%
plot(nsvals,koutvals.*good)
legend({'1','2','3'})
%ylim([0 0.01]);

%% more consistency of fixed kin

kinvals = 10^(-3)*(0.6:0.01:1.3);
nsvals = [];
csvals = [];
koutvals = [];

for k = 1:numel(kinvals)
    kin = kinvals(k)*[1 1 1];
    output{k} = fitKineticModelFixKinFitAlpha(input, kin);
    for i = 1:3
        nsvals(k,i) = output{k}{i}.ns;
        csvals(k,i) = output{k}{i}.cs;
        koutvals(k,i) = output{k}{i}.kout;
    end
end

good = all(koutvals > 0 & (nsvals > 0 & nsvals < 1) & (csvals > 0 & csvals < 1),2);
kinvals(good)

%% visualize model fit 

prefix = ['kinfixed' num2str(num2str(kin(1))) '_'];
outputcomb = cat(2,output{:});
inferred = {'kin','kout','cs','ns'};
vals = {};
errorvals = {};
ylabelstr = {'rate (10^{-3} sec^{-1})', '% change'};
titlestr = {'nuclear exchange','sequestration'};
legendpos = {'NorthEast','NorthWest'};
legendstr = {{'export','import'},{'nuclear','cytoplasm'}};
ylimval = {[0 7],[-100 100]};

% ylabelstr = {'rate (10^{-3} sec^{-1})', 'Smad4 fraction'};
% ylimval = {[0 7],[0 0.6]};

% 
% %%
% 
% j = 2;
% ki = 1;
% 100*(vals{ki}/vals{ki}(1) - 1)
% 100*(vals{ki}/vals{ki}(1)).*sqrt((errorvals{ki}(1)./vals{ki}(1))^2 + (errorvals{ki}./vals{ki}).^2)

for j = 1
    
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
    % this is to get relative changes in sequestration
    % for absolute, comment out
    if j == 2
        for ki = 1:k-1
            errorvals{ki} = 100*(vals{ki}/vals{ki}(1)).*sqrt((errorvals{ki}(1)./vals{ki}(1))^2 + (errorvals{ki}./vals{ki}).^2);
            errorvals{ki}(1) = 0;
%             errorvals{ki} = 100*errorvals{ki}./vals{ki}(1);
            vals{ki} = 100*(vals{ki}./vals{ki}(1) - 1);
        end
    end
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
    saveas(gcf, fullfile(dataDir, fname));
end
