function output = fitKineticModel2(input)
% fit Smad4 kinetic model to FRAP data
% 
% output = fitKineticModel(input)
% 
% input: structure with input parameters

% lower and upper bounds
lb = [0 0 0 0 0];
ub = [1 1 1 1 1];
alpha = 1.2;

kap = @(kin, kout) kin/(kin+kout);
A = @(kin, kout, cs,ns) kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs)/(ns + kap(kin,kout)*(1-ns-cs));
R = @(kin, kout, cs,ns) alpha*(ns + kap(kin,kout)*(1-ns-cs))/(cs + (1-kap(kin,kout))*(1-ns-cs));
k = @(kin, kout) kin/kap(kin,kout);

% weighed vector of differences 
% (E ~ sum f.^2, so terms are weighted by the inverse variance)
f = @(p) [  (A(p(1),p(2),p(3),p(4)) - input.A)/input.sigA,...
            (R(p(1),p(2),p(3),p(4)) - input.R)/input.sigR,...
            (A(p(1),p(5),p(3),p(4)) - input.Ap)/input.sigAp,...
            (R(p(1),p(5),p(3),p(4)) - input.Rp)/input.sigRp,...
            (k(p(1),p(2)) - input.k)/input.sigk,...
            (k(p(1),p(5)) - input.kp)/input.sigkp ];
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

pinit = [kout0 kin0 cs0 ns0 koutp0];

% actual fitting
options.Algorithm = 'trust-region-reflective';
options.TolFun = 1e-10;
options.Display = 'off';
[pfit,~,~,~,~,~,jacobian] = lsqnonlin(f, pinit, lb, ub, options);

output.kin = pfit(1);
output.kout = pfit(2);
output.cs = pfit(3);
output.ns = pfit(4);
output.koutp = pfit(5);

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

output.sigkin = sqrt(Sp(1,1));
output.sigkout = sqrt(Sp(2,2));
output.sigcs = sqrt(Sp(3,3));
output.signs = sqrt(Sp(4,4));
output.sigkoutp = sqrt(Sp(5,5));

disp('inferred parameter values: ');
disp('---------------------------');
disp(['kin:   ' num2str(output.kin,'%.1e') ' (' num2str(output.sigkin,'%.1e') ')']);
disp(['kout:  ' num2str(output.kout,'%.1e') ' (' num2str(output.sigkout,'%.1e') ')']);
disp(['cs:    ' num2str(output.cs,'%.1e') ' (' num2str(output.sigcs,'%.1e') ')']);
disp(['ns:    ' num2str(output.ns,'%.1e') ' (' num2str(output.signs,'%.1e') ')']);
disp(['koutp: ' num2str(output.koutp,'%.1e') ' (' num2str(output.sigkoutp,'%.1e') ')']);

end