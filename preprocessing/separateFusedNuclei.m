function [newNucleiMask, fusedMask] = separateFusedNuclei(nucleiMask, options) 
    % separate fused nuclei in a binary image
    %
    % [newNucleiMask, fusedMask] = separateFusedNuclei(nucleiMask, options) 
    %
    % nucleiMask:       input binary mask of nuclei
    % options:          structure with fields
    % -minAreaStd:      only objects with A > mean(A) + minAreaStd*std(A)
    %                   can be considered fused (default 1)
    % -minSolidity:     only objects with solidity less than this can be
    %                   considered fused (default 0.95)
    %                   NOTE: this part is computationally expensive
    %                   set value <= 0 to turn off and speed up
    %
    % newNucleiMask:    mask with separated nuclei
    % fusedMask:        mask containing potentially fused nuclei
    %
    % uses erosion of fused nuclei as seeds for seeded watershed within the
    % original mask
    
    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------
    
    if ~exist('options','var')
        minSolidity = 0.95;
        minAreaStd = 1;
    else
        minSolidity = options.minSolidity;
        minAreaStd = options.minAreaStd;
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
    
    s = round(sqrt(mean(area))/pi);
    nucmin = imerode(fusedMask,strel('disk',s));
    basin = imcomplement(bwdist(~fusedMask));
    basin = imimposemin(basin, nucmin | ~fusedMask);

    L = watershed(basin);
    newNucleiMask = L > 1 | nucleiMask - fusedMask;
    
    newNucleiMask = bwareaopen(newNucleiMask, round(0.2*mean(area)));
end