function output = fitKineticModel(input)
% fit Smad4 kinetic model to FRAP data
% 
% output = fitKineticModel(input)
% 
% input: structure with input parameters

% lower and upper bounds
lb = [0 0 0 0 0];
ub = [1 1 1 1 1];

A = @(kap,cs,ns) kap*(1-kap)*(1-ns-cs)/(ns + kap*(1-ns-cs));
R = @(kap,cs,ns) (ns + kap*(1-ns-cs))/(cs + (1-kap)*(1-ns-cs));
k = @(kap, kin) kin/kap;

% weighed vector of differences (m means measured)
% (E ~ sum f.^2, so terms are weighted by the inverse variance)
f = @(p) [  (A(p(1),p(3),p(4)) - input.Am)/input.sigAm,...
            (R(p(1),p(3),p(4)) - input.Rm)/input.sigRm,...
            (A(p(2),p(3),p(4)) - input.Apm)/input.sigApm,...
            (R(p(2),p(3),p(4)) - input.Rpm)/input.sigRpm,...
            (k(p(1),p(5)) - input.km)/input.sigkm,...
            (k(p(2),p(5)) - input.kpm)/input.sigkpm ];
% is everything properly normalized? (i.e. it seems to compare different
% parameters the coefficient of variation would be better than the variance
% but also the terms in the energy function should be normalized in that
% case so the normalizations drop out in the weighted sum?)

% initial values
kap0 = 0.5;
kapp0 = 0.5;
cs0 = 0.5;
ns0 = 0.5;
kin0 = 0.001;

pinit = [kap0 kapp0 cs0 ns0 kin0];

% actual fitting
options.Algorithm = 'trust-region-reflective';
options.TolFun = 1e-10;
options.Display = 'off';
[pfit,~,~,~,~,~,jacobian] = lsqnonlin(f, pinit, lb, ub, options);

output.kap = pfit(1);
output.kapp = pfit(2);
output.cs = pfit(3);
output.ns = pfit(4);
output.kin = pfit(5);

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

output.sigkap = sqrt(Sp(1,1));
output.sigkapp = sqrt(Sp(2,2));
output.sigcs = sqrt(Sp(3,3));
output.signs = sqrt(Sp(4,4));
output.sigkin = sqrt(Sp(5,5));

end