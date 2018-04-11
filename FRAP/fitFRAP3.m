function [parameters, frapframe, gof] = fitFRAP3(results)
    % [parameters, frapframe, gof] = fitFRAP(results)
    % 
    % gof = goodness of fit, from function fit
    %
    % function f to be fitted:
    % f = A*(1 - exp(-k t))
    % A recovery fraction
    %
    % or A + B e^(-kt) for bleaching the opposite compartment
    %
    %--------------------------------------------------------
    % v3: A from raw intensities, only fit k
    %--------------------------------------------------------

    if ~isfield(results,'bleachType')
        results.bleachType = 'nuclear';
    end

    if ~isfield(results,'fitType')
        results.fitType = 'nuclear';
    end

    if strcmp(results.fitType,'nuclear')

        tracesnorm = results.tracesNucNorm;
        %tracesnorm = results.tracesNuc' - results.bgempty;
        %tracesnorm = (tracesnorm./max(tracesnorm))';

    elseif strcmp(results.fitType,'cytoplasmic')

        tracesnorm = results.tracesCytNorm;
        %tracesnorm = results.tracesCyt' - results.bgempty;
        %tracesnorm = (tracesnorm./max(tracesnorm))';
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
    Ball = zeros([Nfrapped 3]);
    t = tres*(0:size(tracesnorm,2)-1);
    
    % for tracked nuclei the trace can end before tmax
    tmax = min(tmax, size(tracesnorm,2)); 
    
    % lower and upper bounds
    lb(1) = 0;
    ub(1) = 1;
    lb(2) = 0;
    ub(2) = Inf;
    lb(3) = -1;
    ub(3) = 1;
        
    % when bleached compartment = measured compartment
    % then we are measuring recovery
    if strcmp(results.fitType, results.bleachType)
        disp('fitting recovery');
        func = @(p1,p2,p3,x) p2*(1-exp(-x*p1)) + p3;
        decay = false;
        
    % otherwise decay
    else 
        disp('fitting decay');
        func = @(p1,p2,p3,x) p2*exp(-x*p1) + p3;
        decay = true;
        
    end
    ft = fittype(func, 'problem',{'p2','p3'}); 
    
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
        tdata = t(frapframe:tmaxs) - t(frapframe);
        fdata = tracesnorm(shapeIdx, frapframe:tmaxs);
        
        if ~any(isnan(fdata))

            endval = mean(tracesnorm(shapeIdx, tmaxs-10:tmaxs));
            postfrapval = mean(tracesnorm(shapeIdx, frapframe:frapframe+1));

            if decay
                A0 = postfrapval - endval;
                B0 = endval;
            else
                A0 = endval - postfrapval;
                B0 = postfrapval;
            end
            k0 = 0.01;
            
            % CUTOFF: enforce that the temporal cutoff is within 10% of
            % equilibrium (typically I pick the cutoff when I think it is
            % at equilibrium, but for bad traces the fit can occasionally
            % deviate from that so this is a way to enforce it)
            if isfield(results,'kcutPercent')
                kcut = -log(results.kcutPercent)/(tmaxs*results.tres);
            else
                kcut = -log(0.1)/(tmaxs*results.tres);
            end
            disp(['kcut :' num2str(kcut)]);

            [outfit{shapeIdx}, gof{shapeIdx}] = fit(tdata',fdata',ft,...
                                'Lower',kcut,'Upper',Inf,'StartPoint',k0,...
                                'problem',{A0,B0});

            k = outfit{shapeIdx}.p1;
            A = outfit{shapeIdx}.p2;
            B = outfit{shapeIdx}.p3;

            tau=1/(60*k);

            CI = confint(outfit{shapeIdx})';

            kall(shapeIdx,:) = [k CI(1,:)];
            Aall(shapeIdx,:) = [A 0 0];
            Ball(shapeIdx,:) = [B 0 0];

            parameters = struct('A', Aall, 'B', Ball, 'k', kall);

            disp(['recovery fraction: ' num2str(A,2)]);
            disp(['const: ' num2str(B,2)]);
            disp(['kin + kout: ' num2str(k,2) ' s^-1 --> tau = ' num2str(tau,2) ' min']);
        end
    end

end