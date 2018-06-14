function t = thresholdMP(im)
    % called from processVsi to create binary mask of colonies

    minI = double(min(im(:)));
    maxI = double(max(im(:)));

    im = mat2gray(im);
    t = 0.65*graythresh(im)*maxI + minI; 
end