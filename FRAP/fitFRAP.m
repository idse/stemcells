function [parameters, frapframe, gof] = fitFRAP(results)
    % [parameters, frapframe, gof] = fitFRAP(tracesnorm, frapframe, tmax, tres)
    % 
    % gof = goodness of fit, from function fit
    %
    % function f to be fitted:
    % f = A*(1 - exp(-k t))
    % A recovery fraction
    %
    % or A + B e^(-kt) for bleaching the opposite compartment

    if ~isfield(results,'bleachType')
        results.bleachType = 'nuclear';
    end
    
    if ~isfield(results,'fitType')
        results.fitType = 'nuclear';
    end
    
    if strcmp(results.fitType,'nuclear')
        
        tracesnorm = results.tracesNucNorm;
        
    elseif strcmp(results.fitType,'cytoplasmic')
        
        tracesnorm = results.tracesCytNorm;
    else
        error('unknown fit type');
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
    
    % lower and upper bounds
    lb(1) = 0;
    ub(1) = 1;
    lb(2) = 0;
    ub(2) = Inf;
            
    % when bleached compartment = measured compartment
    % then we are measuring recovery
    if strcmp(results.fitType, results.bleachType)
        disp('fitting recovery');
        func = @(p1,p2,x) p1*(1-exp(-x*p2));
        decay = false;
        
    % otherwise decay
    else 
        disp('fitting decay');
        func = @(p1,p2,p3,x) p1*exp(-x*p2) + p3;
        % bound on extra parameter
        lb(3) = 0;
        ub(3) = 1;
        decay = true;
        
    end
    ft = fittype(func); 
    
    % detect the frap frame if we're looking at the compartment being
    % bleached
    if decay && isfield(results, 'frapframe')
        frapframe = results.frapframe;
    elseif decay && ~isfield(results, 'frapframe')
        warning('frapframe was not defined, using 1');
        frapframe = 1;
    else
        [~,frapframe] = min(tracesnorm(1,1:10));
    end
    
    for shapeIdx = 1:Nfrapped

        if ~isempty(tmax)
            tmaxs = tmax(shapeIdx);
        else
            tmaxs = length(t);
        end
        tdata = t(frapframe:tmaxs);
        fdata = tracesnorm(shapeIdx, frapframe:tmaxs);
        
        if ~any(isnan(fdata))
            
            if decay
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
            if decay
                B = outfit{shapeIdx}.p3;
            end
            tau=1/(60*k);
            
            CI = confint(outfit{shapeIdx})';

            Aall(shapeIdx,:) = [A CI(1,:)];
            kall(shapeIdx,:) = [k CI(2,:)];

            parameters = struct('A', Aall, 'k', kall);
            
            if decay
                Ball(shapeIdx,:) = [B CI(3,:)];
                parameters.B = Ball;
            end
            
            disp(['recovery fraction: ' num2str(A,2)]);
            disp(['kin + kout: ' num2str(k,2) ' s^-1 --> tau = ' num2str(tau,2) ' min']);
        end
    end

end