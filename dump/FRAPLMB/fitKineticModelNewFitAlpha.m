function [output,res] = fitKineticModelNewFitAlpha(input)
% fit Smad4 kinetic model to FRAP data
% 
% new, avoid alpha by not using R
%
% output = fitKineticModel(input)
% 
% input: structure with input parameters

output = {};

% lower and upper bounds
lb = [0 0 0 0 0];
ub = [1 1 1 1 1];

lb = [lb lb lb 0]; 
ub = [ub ub ub Inf];

kap = @(kin, kout) kin/(kin+kout);
An = @(kin, kout, cs, ns) kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs)/(ns + kap(kin,kout)*(1-ns-cs));
Ac = @(kin, kout, cs, ns) kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs)/(cs + (1-kap(kin,kout))*(1-ns-cs));
R = @(kin, kout, cs, ns) (ns + kap(kin,kout)*(1-ns-cs))/(cs + (1-kap(kin,kout))*(1-ns-cs));
k = @(kin, kout) kin/kap(kin,kout);

% weighed vector of differences 
% (E ~ sum f.^2, so terms are weighted by the inverse variance)
f = @(p) [  ( An(p(1),p(2),p(3),p(4))         - input{1}.A)/input{1}.sigA,...
            ( Ac(p(1),p(2),p(3),p(4))      - input{1}.Ac)/input{1}.sigAc,...
            ( k(p(1),p(2))                   - input{1}.k)/input{1}.sigk,...
            ( R(p(1),p(2),p(3),p(4))*p(16)  - input{1}.r)/input{1}.sigr,...
            ...
            ( An(p(1),p(5),p(3),p(4))         - input{1}.Almb)/input{1}.sigAlmb,...
            ( Ac(p(1),p(5),p(3),p(4))      - input{1}.Aclmb)/input{1}.sigAclmb,...
            ( k(p(1),p(5))                   - input{1}.klmb)/input{1}.sigklmb,...
            ( R(p(1),p(5),p(3),p(4))*p(16)  - input{1}.rlmb)/input{1}.sigrlmb,...
            ...
            ( An(p(6),p(7),p(8),p(9))         - input{2}.A)/input{2}.sigA,...
            ( Ac(p(6),p(7),p(8),p(9))      - input{2}.Ac)/input{2}.sigAc,...
            ( k(p(6),p(7))                   - input{2}.k)/input{2}.sigk,...
            ( R(p(6),p(7),p(8),p(9))*p(16)  - input{2}.r)/input{2}.sigr,...
            ...
            ( An(p(6),p(10),p(8),p(9))        - input{2}.Almb)/input{2}.sigAlmb,...
            ( Ac(p(6),p(10),p(8),p(9))    - input{2}.Aclmb)/input{2}.sigAclmb,...
            ( k(p(6),p(10))                  - input{2}.klmb)/input{2}.sigklmb,...
            ( R(p(6),p(10),p(8),p(9))*p(16)  - input{2}.rlmb)/input{2}.sigrlmb,...
            ...
            ( An(p(11),p(12),p(13),p(14))     - input{3}.A)/input{3}.sigA,...
            ( Ac(p(11),p(12),p(13),p(14))  - input{3}.Ac)/input{3}.sigAc,...
            ( k(p(11),p(12))                 - input{3}.k)/input{3}.sigk,...
            ( R(p(11),p(12),p(13),p(14))*p(16)  - input{3}.r)/input{2}.sigr,...
            ...
            ( An(p(11),p(15),p(13),p(14))     - input{3}.Almb)/input{3}.sigAlmb,...
            ( Ac(p(11),p(15),p(13),p(14))  - input{3}.Aclmb)/input{3}.sigAclmb,...
            ( k(p(11),p(15))                 - input{3}.klmb)/input{3}.sigklmb,...
            ( R(p(11),p(15),p(13),p(14))*p(16)  - input{3}.rlmb)/input{2}.sigrlmb];

errorvec = [];
for i = 1:3
    errorvec = [errorvec input{i}.sigA input{i}.sigAc input{i}.sigk input{i}.sigr...
                input{i}.sigAlmb input{i}.sigAclmb input{i}.sigklmb input{i}.sigrlmb];
end

% is everything properly normalized? (i.e. it seems to compare different
% parameters the coefficient of variation would be better than the variance
% but also the terms in the energy function should be normalized in that
% case so the normalizations drop out in the weighted sum?)

% initial values
kout0 = 0.003;
koutp0 = kout0/10;
cs0 = 0.5;
ns0 = 0.5;
kin0 = 0.001;
alpha0 = 1;

pinit = [kout0 kin0 cs0 ns0 koutp0];
pinit = [pinit pinit pinit alpha0]; 
 
% actual fitting
options.Algorithm = 'trust-region-reflective';
options.TolFun = 1e-10;
options.Display = 'off';
[pfit,resnorm,res,~,~,~,jacobian] = lsqnonlin(f, pinit, lb, ub, options);

for i = 1:3
    output{i}.kin = pfit(5*(i-1)+1);
    output{i}.kout = pfit(5*(i-1)+2);
    output{i}.cs = pfit(5*(i-1)+3);
    output{i}.ns = pfit(5*(i-1)+4);
    output{i}.koutp = pfit(5*(i-1)+5);
    output{i}.alpha = pfit(16);
end

%error('test');

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
%actualJ = repmat(errorvec', [1 size(J,2)]).*J;

disp('---------fit +alpha, fit to A_c and R------------');
disp(['residual norm: ' num2str(resnorm) ' sqrt = ' num2str(sqrt(resnorm))]);

for i = 1:3

    k = 5*(i-1);
    output{i}.sigkin = sqrt(Sp(k+1,k+1));
    output{i}.sigkout = sqrt(Sp(k+2,k+2));
    output{i}.sigcs = sqrt(Sp(k+3,k+3));
    output{i}.signs = sqrt(Sp(k+4,k+4));
    output{i}.sigkoutp = sqrt(Sp(k+5,k+5));

    disp('inferred parameter values: ');
    disp('---------------------------');
    disp(['kin:   ' num2str(output{i}.kin,'%.1e') ' (' num2str(output{i}.sigkin,'%.1e') ')']);
    disp(['kout:  ' num2str(output{i}.kout,'%.1e') ' (' num2str(output{i}.sigkout,'%.1e') ')']);
    disp(['cs:    ' num2str(output{i}.cs,'%.1e') ' (' num2str(output{i}.sigcs,'%.1e') ')']);
    disp(['ns:    ' num2str(output{i}.ns,'%.1e') ' (' num2str(output{i}.signs,'%.1e') ')']);
    disp(['koutp: ' num2str(output{i}.koutp,'%.1e') ' (' num2str(output{i}.sigkoutp,'%.1e') ')']);
end
disp(['alpha:    ' num2str(output{1}.alpha,'%.1e')]);

end