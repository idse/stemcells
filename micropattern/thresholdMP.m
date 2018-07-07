function t = thresholdMP(im, adjustmentFactor)
    % called from processVsi to create binary mask of colonies

    if ~exist('adjustmentFactor','var') || isempty(adjustmentFactor)
        adjustmentFactor = 0.65;
    end
    
    minI = double(min(im(:)));
    maxI = double(max(im(:)));

    im = mat2gray(im);
    t = adjustmentFactor*graythresh(im)*maxI + minI; 
end