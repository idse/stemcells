%% manual checks

Ninit = {[191.8	231.0	101.6; 120.8	244.3	107.4]};
Cinit = {[167.7	454.1	176.4; 128.2	539.8	198.4]};

Nbleach = {[204.0	144.1	47.9; 106.9	139.5	43.7]};
Cbleach = {[166.8	363.8	148.0; 125.0	422.9	166.2]};

Nfin = {[219.1	167.4	68.5; 98.9	162.4	64.4]};
Cfin = {[114.9	328.6	137.9; 96.3	363.6	151.1]};

bg = {[301.4	146.5	49.3]};
ctrlinit = {[10094.0	356.5	148.9]};
ctrlfin = {[7828.0	318.6	137.2]};

BF = 0.83;

i = 1;
vidx = [2 1];

disp('----')

bgval = min(cat(1, bg{i}(:,2), Nbleach{i}(:,2)));
BF = mean((ctrlfin{i}(:,2) - bgval)./(ctrlinit{i}(:,2) - bgval));

R = (Ninit{i}(:,2) - bgval)./(Cinit{i}(:,2) - bgval);
Rp = (Nfin{i}(:,2) - bgval)./(Cfin{i}(:,2) - bgval); % no BF there: it drops out
Rrat = Rp./R;

%disp('Rrat')
%disp([Rrat allOut{1}.Rrat(vidx)]);

nb = (Nbleach{i}(:,2) - bgval)./(Ninit{i}(:,2) - bgval);
nc = (Cbleach{i}(:,2) - bgval)./(Cinit{i}(:,2) - bgval);
corrfac = nc - nb;

Araw = ((Nfin{i}(:,2)-bgval)/BF + bgval - Nbleach{i}(:,2))./(Ninit{i}(:,2) - bgval);
ArawCyt = (Cbleach{i}(:,2) - (Cfin{i}(:,2)-bgval)/BF -bgval)./(Cinit{i}(:,2) - bgval);

Acorr = Araw./corrfac; 
%disp('Acorr');
%disp([Acorr allOut{1}.Acorr(vidx)]);

AcorrCyt = ArawCyt./corrfac; 
%disp('AcorrCyt');
%disp([AcorrCyt allOut{1}.AcorrCyt(vidx)]);

%Ninitval = [Ninit{i}(:,2) allOut{1}.Ninit(vidx)];
%Nfinval = [Nfin{i}(:,2) allOut{1}.Nfin(vidx)];
%Cinitval = [Cinit{i}(:,2) allOut{1}.Cinit(vidx)];
%Cfinval = [Cfin{i}(:,2) allOut{1}.Cfin(vidx)];
