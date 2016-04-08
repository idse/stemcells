function [radialMaskStack, edges] = makeRadialBinningMasks(meta)
    % returns a stack of radial masks for radial averaging w/o segmenting
    %
    % [radialMaskStack, edges] = makeRadialBinningMasks(meta)
    %
    % meta:     metadata object
    % edges:    cell array of bin edges (5 micron apart) for different size
    %           colonies
    % radialMaskStack:  cell array of stacks of radial mask
    
    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------

    radialMaskStack ={};
    edges = {};
    
    Rpixel = ceil(meta.colRadiiPixel);
    
    % I want 5 micron wide bins (~1/2 cell diameter)
    binWidthMicron = 5;
    
    for i = 1:numel(Rpixel)

        Rcol = Rpixel(i);
        Rmax = Rcol + meta.colMargin;

        N = meta.colRadiiMicron(i)/binWidthMicron;
        
        [x,y] = meshgrid(-Rmax:Rmax,-Rmax:Rmax);
        r = sqrt(x.^2 + y.^2);

        edges{i} = sqrt(linspace(0,Rmax^2,N+1));
        radialMaskStack{i} = false([size(x) N]);
        for ri = 1:N
            radialMaskStack{i}(:,:,ri) = r > edges{i}(ri) & r < edges{i}(ri+1);
        end
    end 
end