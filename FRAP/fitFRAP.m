function [Aall, kall, frapframe] = fitFRAP(tracesnorm, tmax, tres)
    % [Aall, kall] = fitFRAP(tracesnorm, frapframe, tmax, tres)
    % 
    % function f to be fitted:
    % f = A*(1 - exp(-k t))
    % A recovery fraction

    outfit = {};
    Nfrapped = size(tracesnorm,1);
    Aall = zeros([Nfrapped 3]);
    kall = zeros([Nfrapped 3]);
    t = tres*(0:size(tracesnorm,2)-1);
    
    % detect the frap frame
    [~,frapframe] = min(tracesnorm(1,:));
    
    func = @(p1,p2,x) p1*(1-exp(-x*p2));
    ft = fittype(func); 
        
    for shapeIdx = 1:Nfrapped

        tdata = t(frapframe:tmax(shapeIdx));
        fdata = tracesnorm(shapeIdx, frapframe:tmax(shapeIdx));

        R0 = 0.5;
        k0 = 0.01;

        pinit = [R0 k0];

        % lower and upper bounds
        lb(1) = 0;
        ub(1) = 1;
        lb(2) = 0;
        ub(2) = Inf;

        outfit{shapeIdx} = fit(tdata',fdata',ft,'Lower',lb,'Upper',ub,'StartPoint',pinit);

        A = outfit{shapeIdx}.p1;
        k = outfit{shapeIdx}.p2;
        tau=1/(60*k);

        CI = confint(outfit{shapeIdx})';

        Aall(shapeIdx,:) = [A CI(1,:)];
        kall(shapeIdx,:) = [k CI(2,:)];

        disp(['recovery fraction: ' num2str(A,2)]);
        disp(['kin + kout: ' num2str(k,2) ' s^-1 --> tau = ' num2str(tau,2) ' min']);
    end

end