function [Aall, kall, frapframe, gof] = fitFRAP(results)
    % [Aall, kall, frapframe, gof] = fitFRAP(tracesnorm, frapframe, tmax, tres)
    % 
    % gof = goodness of fit, from function fit
    %
    % function f to be fitted:
    % f = A*(1 - exp(-k t))
    % A recovery fraction
    %
    % or A + B e^(-kt) for bleaching the opposite compartment

    if ~isfield(results,'bleachType')
        bleachType = 'nuclear';
    else
        bleachType = results.bleachType;
    end
    
    if ~isfield(results,'fitType')
        results.fitType = 'recovery';
    end
    
    if strcmp(bleachType,'nuclear')
        tracesnorm = results.tracesNucNorm;
    elseif strcmp(bleachType,'cytoplasmic')
        tracesnorm = results.tracesCytNorm;
    else
        error('unknown bleach type');
    end

    tmax = results.tmax;
    tres = results.tres;
    
    outfit = {};
    gof = {};
    Nfrapped = size(tracesnorm,1);
    Aall = zeros([Nfrapped 3]);
    kall = zeros([Nfrapped 3]);
    t = tres*(0:size(tracesnorm,2)-1);
    
    % for tracked nuclei the trace can end before tmax
    tmax = min(tmax, size(tracesnorm,2)); 
    
    % detect the frap frame
    [~,frapframe] = min(tracesnorm(1,1:10));
    
    % lower and upper bounds
    lb(1) = 0;
    ub(1) = 1;
    lb(2) = 0;
    ub(2) = Inf;
            
    if strcmp(results.fitType, 'decay')
        func = @(p1,p2,p3,x) p1*exp(-x*p2) + p3;
        % bound on extra parameter
        lb(3) = 0;
        ub(3) = 1;
    else 
        func = @(p1,p2,x) p1*(1-exp(-x*p2));
    end
    ft = fittype(func); 
    
    for shapeIdx = 1:Nfrapped

        tdata = t(frapframe:tmax(shapeIdx));
        fdata = tracesnorm(shapeIdx, frapframe:tmax(shapeIdx));
        
        if ~any(isnan(fdata))
            
            if strcmp(results.fitType, 'decay')
                A0 = 1;
                k0 = 0.01;
                B0 = 0;
                pinit = [A0 k0 B0];
            else
                A0 = 0.5;
                k0 = 0.01;
                pinit = [A0 k0];
            end

            [outfit{shapeIdx}, gof{shapeIdx}] = fit(tdata',fdata',ft,'Lower',lb,'Upper',ub,'StartPoint',pinit);

            A = outfit{shapeIdx}.p1;
            k = outfit{shapeIdx}.p2;
            if strcmp(results.fitType, 'decay')
                B = outfit{shapeIdx}.p3;
            end
            tau=1/(60*k);

            CI = confint(outfit{shapeIdx})';

            Aall(shapeIdx,:) = [A CI(1,:)];
            kall(shapeIdx,:) = [k CI(2,:)];
            
            disp(['recovery fraction: ' num2str(A,2)]);
            disp(['kin + kout: ' num2str(k,2) ' s^-1 --> tau = ' num2str(tau,2) ' min']);
        end
    end

end