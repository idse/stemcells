function [newNucleiMask, fusedMask] = separateFusedNuclei(nucleiMask, options) 
    % separate fused nuclei in a binary image
    %
    % [newNucleiMask, fusedMask] = separateFusedNuclei(nucleiMask, options) 
    %
    % uses erosion of fused nuclei as seeds for seeded watershed within the
    % original mask
    %
    % nucleiMask:       input binary mask of nuclei
    % options:          structure with fields
    % -minAreaStd:      only objects with A > mean(A) + minAreaStd*std(A)
    %                   can be considered fused (default 1)
    % -minSolidity:     only objects with solidity less than this can be
    %                   considered fused (default 0.95)
    %                   NOTE: this part is computationally expensive
    %                   set value <= 0 to turn off and speed up
    % -erodeSize        in units of mean radius, default 1
    %
    % newNucleiMask:    mask with separated nuclei
    % fusedMask:        mask containing potentially fused nuclei

    
    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------
    
    if ~exist('options','var')
        options = struct();
    end
    if ~isfield(options,'minSolidity')
        minSolidity = 0.95;
    else
        minSolidity = options.minSolidity;
    end
    if ~isfield(options,'minAreaStd')
        minAreaStd = 1;
    else
        minAreaStd = options.minAreaStd;
    end
    if ~isfield(options,'erodeSize')
        erodeSize = 1;
    else
        erodeSize = options.erodeSize;
    end
    
    CC = bwconncomp(nucleiMask);
    if minSolidity > 0
        stats = regionprops(CC, 'ConvexArea', 'Area');
        convexArea = [stats.ConvexArea];
    else
        stats = regionprops(CC, 'Area');
    end
    area = [stats.Area];
    
    if minSolidity > 0
        fusedCandidates = area./convexArea < minSolidity & area > mean(area) + minAreaStd*std(area);
    else
        fusedCandidates = area > mean(area) + minAreaStd*std(area);
    end
    sublist = CC.PixelIdxList(fusedCandidates);
    sublist = cat(1,sublist{:});

    fusedMask = false(size(nucleiMask));
    fusedMask(sublist) = 1;

%     figure,
%     imshow(cat(3,nucleiMask,fused,0*fused))

    s = round(erodeSize*sqrt(mean(area))/pi);
    nucmin = imerode(fusedMask,strel('disk',s));
    outside = ~imdilate(fusedMask,strel('disk',1));
    basin = imcomplement(bwdist(outside));
    basin = imimposemin(basin, nucmin | outside);

    L = watershed(basin);
    newNucleiMask = L > 1 | nucleiMask - fusedMask;
    
    newNucleiMask = bwareaopen(newNucleiMask, round(0.2*mean(area)));
end