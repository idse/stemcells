function output = fitKineticModel2(input, alpha)
% fit Smad4 kinetic model to FRAP data
% 
% output = fitKineticModel(input)
% 
% input: structure with input parameters

output = {};

% lower and upper bounds
lb = [0 0 0 0 0];
ub = [1 1 1 1 1];

for i = 1:3
    
    inputi = input{i};

    kap = @(kin, kout) kin/(kin+kout);
    A = @(kin, kout, cs, ns) kap(kin,kout)*(1-kap(kin,kout))*(1-ns-cs)/(ns + kap(kin,kout)*(1-ns-cs));
    R = @(kin, kout, cs, ns) alpha*(ns + kap(kin,kout)*(1-ns-cs))/(cs + (1-kap(kin,kout))*(1-ns-cs));
    k = @(kin, kout) kin/kap(kin,kout);

    % weighed vector of differences 
    % (E ~ sum f.^2, so terms are weighted by the inverse variance)
    f = @(p) [  (A(p(1),p(2),p(3),p(4)) - inputi.A)/inputi.sigA,...
                (R(p(1),p(2),p(3),p(4)) - inputi.R)/inputi.sigR,...
                (A(p(1),p(5),p(3),p(4)) - inputi.Ap)/inputi.sigAp,...
                (R(p(1),p(5),p(3),p(4)) - inputi.Rp)/inputi.sigRp,...
                (k(p(1),p(2)) - inputi.k)/inputi.sigk,...
                (k(p(1),p(5)) - inputi.kp)/inputi.sigkp ];
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
    %alpha = 1.2;

    pinit = [kout0 kin0 cs0 ns0 koutp0];

    % actual fitting
    options.Algorithm = 'trust-region-reflective';
    options.TolFun = 1e-10;
    options.Display = 'off';
    [pfit,~,~,~,~,~,jacobian] = lsqnonlin(f, pinit, lb, ub, options);

    outputi = struct();
    outputi.kin = pfit(1);
    outputi.kout = pfit(2);
    outputi.cs = pfit(3);
    outputi.ns = pfit(4);
    outputi.koutp = pfit(5);
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

    outputi.sigkin = sqrt(Sp(1,1));
    outputi.sigkout = sqrt(Sp(2,2));
    outputi.sigcs = sqrt(Sp(3,3));
    outputi.signs = sqrt(Sp(4,4));
    outputi.sigkoutp = sqrt(Sp(5,5));
    %output.sigalpha = sqrt(Sp(6,6));

    output{i} = outputi;
    
    disp('inferred parameter values: ');
    disp('---------------------------');
    disp(['kin:   ' num2str(outputi.kin,'%.1e') ' (' num2str(outputi.sigkin,'%.1e') ')']);
    disp(['kout:  ' num2str(outputi.kout,'%.1e') ' (' num2str(outputi.sigkout,'%.1e') ')']);
    disp(['cs:    ' num2str(outputi.cs,'%.1e') ' (' num2str(outputi.sigcs,'%.1e') ')']);
    disp(['ns:    ' num2str(outputi.ns,'%.1e') ' (' num2str(outputi.signs,'%.1e') ')']);
    disp(['koutp: ' num2str(outputi.koutp,'%.1e') ' (' num2str(outputi.sigkoutp,'%.1e') ')']);
    %disp(['alpha: ' num2str(output.alpha,'%.1e') ' (' num2str(output.sigalpha,'%.1e') ')']);
end

end