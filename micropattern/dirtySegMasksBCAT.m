function [nonmembrane, membrane] = dirtySegMasksBCAT(colonymask, mask, param)
    % [nuclear, cytoplasmic] = dirtySegMasksBCAT(cleanmask, mask, param)
    % param: struct with parameters
    %   cytoSize, nucShrinkage, minArea
    
    mask = imdilate(mask,strel('disk',1));
end