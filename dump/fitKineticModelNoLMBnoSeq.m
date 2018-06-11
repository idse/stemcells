function output = fitKineticModelNoLMBnoSeq(input)
% fit Smad4 kinetic model to FRAP data
% 
% new, avoid alpha by not using R
% this version, hold some parameter fixed to avoid LMB data
%
% output = fitKineticModel(input)
% 
% input: structure with input parameters

output = {};

% lower and upper bounds
lb = [-1 -1];
ub = Inf*[1 1];

lb = [lb lb lb 1]; 
ub = [ub ub ub Inf];
    
kap = @(kin, kout) kin/(kin+kout);
An = @(kin, kout, cs, ns) kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs)/(ns + kap(kin,kout)*(1-ns-cs));
Ac = @(kin, kout, cs, ns) kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs)/(cs + (1-kap(kin,kout))*(1-ns-cs));
R = @(kin, kout, cs, ns) (ns + kap(kin,kout)*(1-ns-cs))/(cs + (1-kap(kin,kout))*(1-ns-cs));
k = @(kin, kout) kin/kap(kin,kout);

ns = [0 0 0];
cs = [0 0 0];

% weighed vector of differences 
% (E ~ sum f.^2, so terms are weighted by the inverse variance)
f = @(p) [  ( An(p(1),p(2),cs(1),ns(1))      - input{1}.A)/input{1}.sigA,...
            ( Ac(p(1),p(2),cs(1),ns(1))      - input{1}.Ac)/input{1}.sigAc,...
            ( k(p(1),p(2))                  - input{1}.k)/input{1}.sigk,...
            ( R(p(1),p(2),cs(1),ns(1))*p(7) - input{1}.R)/input{1}.sigR,...
            ( An(p(3),p(4),cs(2),ns(2))      - input{2}.A)/input{2}.sigA,...
            ( Ac(p(3),p(4),cs(2),ns(2))      - input{2}.Ac)/input{2}.sigAc,...
            ( k(p(3),p(4))                  - input{2}.k)/input{2}.sigk,...
            ( R(p(3),p(4),cs(2),ns(2))*p(7) - input{2}.R)/input{2}.sigR,...
            ( An(p(5),p(6),cs(3),ns(3))      - input{3}.A)/input{3}.sigA,...
            ( Ac(p(5),p(6),cs(3),ns(3))      - input{3}.Ac)/input{3}.sigAc,...
            ( k(p(5),p(6))                  - input{3}.k)/input{3}.sigk,...
            ( R(p(5),p(6),cs(3),ns(3))*p(7) - input{3}.R)/input{3}.sigR];

% is everything properly normalized? (i.e. it seems to compare different
% parameters the coefficient of variation would be better than the variance
% but also the terms in the energy function should be normalized in that
% case so the normalizations drop out in the weighted sum?)

% initial values
kout0 = 0.003;
kin0 = 0.001;
alpha0 = 1;

pinit = [kout0 kin0];
pinit = [pinit pinit pinit alpha0]; 

% actual fitting
options.Algorithm = 'trust-region-reflective';
options.TolFun = 1e-10;
options.Display = 'off';
[pfit,resnorm,~,~,~,~,jacobian] = lsqnonlin(f, pinit, lb, ub, options);

% % alternative
% ft = fittype(f, 'problem',{'p2','p3'}); 
% [outfit, gof] = fit(tdata',fdata',ft,...
%                                 'Lower',kcut,'Upper',Inf,'StartPoint',k0,...
%                                 'problem',{A0,B0});

for i = 1:3
    output{i}.kin = pfit(2*(i-1)+1);
    output{i}.kout = pfit(2*(i-1)+2);
    output{i}.cs = cs(i);
    output{i}.ns = ns(i);
    output{i}.alpha = pfit(7);
end

%error('test');
%output.alpha = pfit(6);

% error bars 
%
% E~ sum_i f_i^2 = sum_i dx^T S^(-1) dx = sum_i (A_i - Am_i)^2/sigAm_i^2 
% where dx is the difference vector (dA, dR, dAp, dRp)
% and the covariance matrix S = diag(sigA^2, sigR^2, sigAp^2, sigRp^2)
%
% since S = <dx dx>, a coordinate transformation takes it to 
% Sp = <dxp dxp> = <J dx J dx> = J S J^T
% where J is the Jacobian: J_ij = \partial xp_i / \partial x_j
% and xp is the vector of new parameters (kap, kapp, cs, ns)
%
% HOWEVER, the Jacobian spit out be lsqnonlin is the inverse
% J = \partial f_i / \partial xp_j = \partial x_i/\partial xp_j / sigx_i
% so here dx_i = sigx_i J_ij dxp_j , and
% dx^T S^(-1) dx -> dxp^T J^T J dxp, 
% where the sigx_i^2 dropped out because S^-1 is a diagonal matrix with
% 1/sigx_i^2 on the diagonal.
% 
% that means the diagonal elements of (J^T J)^-1 are the new variances

J = full(jacobian);
Sp = inv(J'*J);
disp('---------fitKineticModel No LMB, No sequestration------------');
disp(['residual norm: ' num2str(resnorm)]);

for i = 1:3

    k = 2*(i-1);
    output{i}.sigkin = sqrt(Sp(k+1,k+1));
    output{i}.sigkout = sqrt(Sp(k+2,k+2));
    %output{i}.sigcs = 0;%sqrt(Sp(k+3,k+3));
    %output{i}.signs = sqrt(Sp(k+3,k+3));
%    output{i}.sigkoutp = sqrt(Sp(k+5,k+5));

    disp('inferred parameter values: ');
    disp('---------------------------');
    disp(['kin:   ' num2str(output{i}.kin,'%.1e') ' (' num2str(output{i}.sigkin,'%.1e') ')']);
    disp(['kout:  ' num2str(output{i}.kout,'%.1e') ' (' num2str(output{i}.sigkout,'%.1e') ')']);
    %disp(['cs:    ' num2str(output{i}.cs,'%.1e') ' (' num2str(output{i}.sigcs,'%.1e') ')']);
    %disp(['ns:    ' num2str(output{i}.ns,'%.1e') ' (' num2str(output{i}.signs,'%.1e') ')']);
%    disp(['koutp: ' num2str(output{i}.koutp,'%.1e') ' (' num2str(output{i}.sigkoutp,'%.1e') ')']);
end
disp(['alpha:    ' num2str(output{1}.alpha,'%.1e')]);
end