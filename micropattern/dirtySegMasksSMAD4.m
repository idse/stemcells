function [nuclear, cytoplasmic] = dirtySegMasksSMAD4(colonymask, mask, param)
    % [nuclear, cytoplasmic] = dirtySegMasksSMAD4(cleanmask, mask, param)
    % param: struct with parameters
    %   cytoSize, nucShrinkage, minArea

    nuclear = colonymask & mask;    
    cytoplasmic = colonymask & (imdilate(nuclear,strel('disk',param.cytoSize))...
                                    - nuclear); 
    nuclear = imerode(nuclear,strel('disk',param.nucShrinkage));
    nuclear = bwareaopen(nuclear,param.minArea);
end