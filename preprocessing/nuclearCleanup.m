function newNuclearMask = nuclearCleanup(nuclearMask, options)
    % clean up nuclear mask
    %
    % newNuclearMask = nuclearCleanup(nuclearMask, options)
    %
    % nuclearMask:      input binary mask of nuclei
    %
    % newNucleiMask:    clean mask
    %    
    % options:          structure with fields:
    %
    % -minArea          remove junk smaller than this (in pixels)
    % -openSize         disk radius for imopen strel
    %
    % -separateFused    boolean
    % -clearBorder      boolean
    %
    % options for separateFusedNuclei:
    %
    % -minAreaStd:      only objects with A > mean(A) + minAreaStd*std(A)
    %                   can be considered fused (default 1)
    % -minSolidity:     only objects with solidity less than this can be
    %                   considered fused (default 0.95)
    %                   NOTE: this part is computationally expensive
    %                   set value <= 0 to turn off and speed up
    % -erodeSize        in units of mean radius, default 1
    
    
    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------

    if ~exist('options','var')
        options = struct();
    end
    if ~isfield(options,'minArea')
        options.minArea = 500;
    end
    if ~isfield(options,'openSize')
        options.openSize = 5;
    end
    if ~isfield(options,'separateFused')
        options.separateFused = true;
    end
    if ~isfield(options,'clearBorder')
        options.clearBorder = true;
    end
    
    nuclearMask = bwareaopen(nuclearMask, options.minArea);    
    
    nuclearMask = imopen(nuclearMask, strel('disk',options.openSize));
    nuclearMask = imfill(nuclearMask,'holes');
    
    if options.separateFused
        nuclearMask = separateFusedNuclei(nuclearMask,options);
    end
    if options.clearBorder 
        newNuclearMask = imclearborder(nuclearMask);
    end
end