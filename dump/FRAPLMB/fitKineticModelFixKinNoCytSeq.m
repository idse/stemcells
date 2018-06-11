function output = fitKineticModelFixKinNoCytSeq(input)
% fit Smad4 kinetic model to FRAP data
% 
% version 3: fit all 3 conditions at once and fit a single alpha
%
% output = fitKineticModel(input)
% 
% input: structure with input parameters

output = {};

% lower and upper bounds
lb = [0 0 0];
ub = [1 1 1];

lb = [0 lb lb lb 0.5]; % first is kin, last value is for alpha
ub = [1 ub ub ub 2];
    
kap = @(kin, kout) kin/(kin+kout);
A = @(kin, kout, ns) kap(kin,kout)*(1-kap(kin,kout))*(1-ns)/(ns + kap(kin,kout)*(1-ns));
R = @(kin, kout, ns, alpha) alpha*(ns + kap(kin,kout)*(1-ns))/(1-kap(kin,kout))*(1-ns);
k = @(kin, kout) kin/kap(kin,kout);

W = [];
for i = 1:3
    W = [W  1/input{i}.sigA 1/input{i}.sigR 1/input{i}.sigAp...
            1/input{i}.sigRp 1000/input{i}.sigk 1000/input{i}.sigkp];
end
   
% weighed vector of differences 
% (E ~ sum f.^2, so terms are weighted by the inverse variance)
f = @(p) [  W(1)*(A(p(1),p(2),p(3))         - input{1}.A),...
            W(2)*(R(p(1),p(2),p(3),p(11))   - input{1}.R),...
            W(3)*(A(p(1),p(4),p(3))         - input{1}.Ap),...
            W(4)*(R(p(1),p(4),p(3),p(11))   - input{1}.Rp),...
            W(5)*(k(p(1),p(2))            - input{1}.k),...
            W(6)*(k(p(1),p(5))            - input{1}.kp),...
            W(7)*(A(p(1),p(5),p(6))         - input{2}.A),...
            W(8)*(R(p(1),p(5),p(6),p(11))   - input{2}.R),...
            W(9)*(A(p(1),p(7),p(6))         - input{2}.Ap),...
            W(10)*(R(p(1),p(7),p(6),p(11))   - input{2}.Rp),...
            W(11)*(k(p(1),p(5))            - input{2}.k),...
            W(12)*(k(p(1),p(7))            - input{2}.kp),...
            W(13)*(A(p(1),p(8),p(9))         - input{3}.A),...
            W(14)*(R(p(1),p(8),p(9),p(11))   - input{3}.R),...
            W(15)*(A(p(1),p(10),p(9))        - input{3}.Ap),...
            W(16)*(R(p(1),p(10),p(9),p(11))  - input{3}.Rp),...
            W(17)*(k(p(1),p(8))                 - input{3}.k),...
            W(18)*(k(p(1),p(10))                 - input{3}.kp)];

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

pinit = [kout0 ns0 koutp0];
pinit = [kin0 pinit pinit pinit alpha0]; 
 
% actual fitting
options.Algorithm = 'trust-region-reflective';
options.TolFun = 1e-10;
options.Display = 'off';
[pfit,~,~,~,~,~,jacobian] = lsqnonlin(f, pinit, lb, ub, options);

for i = 1:3
    output{i}.kin = pfit(1);
    output{i}.kout = pfit(3*(i-1)+2);
    output{i}.cs = 0;%pfit(3*(i-1)+3);
    output{i}.ns = pfit(3*(i-1)+3);
    output{i}.koutp = pfit(3*(i-1)+4);
    output{i}.alpha = pfit(11);
end

W = ones([18 1]);
f = @(p) [  W(1)*(A(p(1),p(2),p(3))          - input{1}.A),...
            W(2)*(R(p(1),p(2),p(3),p(11))    - input{1}.R),...
            W(3)*(A(p(1),p(4),p(3))          - input{1}.Ap),...
            W(4)*(R(p(1),p(4),p(3),p(11))    - input{1}.Rp),...
            W(5)*(k(p(1),p(2))               - input{1}.k),...
            W(6)*(k(p(1),p(5))               - input{1}.kp),...
            W(7)*(A(p(1),p(5),p(6))          - input{2}.A),...
            W(8)*(R(p(1),p(5),p(6),p(11))    - input{2}.R),...
            W(9)*(A(p(1),p(7),p(6))          - input{2}.Ap),...
            W(10)*(R(p(1),p(7),p(6),p(11))   - input{2}.Rp),...
            W(11)*(k(p(1),p(5))              - input{2}.k),...
            W(12)*(k(p(1),p(7))              - input{2}.kp),...
            W(13)*(A(p(1),p(8),p(9))         - input{3}.A),...
            W(14)*(R(p(1),p(8),p(9),p(11))   - input{3}.R),...
            W(15)*(A(p(1),p(10),p(9))        - input{3}.Ap),...
            W(16)*(R(p(1),p(10),p(9),p(11))  - input{3}.Rp),...
            W(17)*(k(p(1),p(8))              - input{3}.k),...
            W(18)*(k(p(1),p(10))             - input{3}.kp)];
%%f(pfit)

res = sqrt(f(pfit).^2);

for i = 1:3
    output{i}.resA = res(6*(i-1)+1);
    output{i}.resR = res(6*(i-1)+2);
    output{i}.resAp = res(6*(i-1)+3);
    output{i}.resRp = res(6*(i-1)+4);
    output{i}.resk = res(6*(i-1)+5);
    output{i}.reskp = res(6*(i-1)+6);
    output{i}.restot = sum(res);
end

%R(pfit(1),pfit(7),pfit(6),pfit(11))
%error('test')

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
%error('test');
disp('---------fitKineticModelFixKinNoCytSeq------------');
for i = 1:3

    k = 1 + 3*(i-1);
    output{i}.sigkin = sqrt(Sp(1,1));
    output{i}.sigkout = sqrt(Sp(k+1,k+1));
    output{i}.sigcs = 0;% sqrt(Sp(k+2,k+2));
    output{i}.signs = sqrt(Sp(k+2,k+2));
    output{i}.sigkoutp = sqrt(Sp(k+3,k+3));
    output{i}.sigalpha = sqrt(Sp(11,11));

    disp('inferred parameter values: ');
    disp('---------------------------');
    disp(['kin:   ' num2str(output{i}.kin,'%.1e') ' (' num2str(output{i}.sigkin,'%.1e') ')']);
    disp(['kout:  ' num2str(output{i}.kout,'%.1e') ' (' num2str(output{i}.sigkout,'%.1e') ')']);
    %disp(['cs:    ' num2str(output{i}.cs,'%.1e') ' (' num2str(output{i}.sigcs,'%.1e') ')']);
    disp(['ns:    ' num2str(output{i}.ns,'%.1e') ' (' num2str(output{i}.signs,'%.1e') ')']);
    disp(['koutp: ' num2str(output{i}.koutp,'%.1e') ' (' num2str(output{i}.sigkoutp,'%.1e') ')']);
    %disp(['alpha: ' num2str(output.alpha,'%.1e') ' (' num2str(output.sigalpha,'%.1e') ')']);
end
disp(['alpha: ' num2str(output{i}.alpha,'%.1e') ' (' num2str(output{i}.sigalpha,'%.1e') ')']);

end