function [divIdx, divAxis] = decideDivide(g, dt, rule, param)
    % decide which cells will divide
    %
    % [divIdx, divAxis] = decideDivide(g, tau, rule, param)
    %
    % divIdx:   binary index of dividing cells
    % divAxis:  division axis
    %
    % g:        lattice
    % dt:      unit of time (rate of division)
    % rule:     name of division law
    % param:    rule parameters
    %
    % 'random', 'area'
    
    % pick a random cell to divide
    %pdivide = ones([1 Ncells])/Ncells;
    %ix = 1:Ncells;

    % nearest neighbor number
    nnn = cellfun(@numel, g.cells);

    S = [numel(g.cells) 1];
    
    % division probability per unit time 
    if      strcmp(rule, 'random')
        
        p = 1/2*ones(S);
        
    elseif  strcmp(rule, 'area')
       
        A0 = param.A0;
        k = param.k;
        a = g.A/A0;
        p = ((nnn' > 4).*(a > 1).*(a-1).^2)./(k+(a-1).^2);
    end

    % probability given the time interval
    P = 1 - (1-p).^dt;

    % division decision
    divIdx = rand(S) < P(:);

    % pick a random division axis
    phi = rand(S)*2*pi;
    divAxis = [cos(phi) sin(phi)];
end