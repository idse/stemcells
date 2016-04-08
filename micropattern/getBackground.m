function bg = getBackground(bg, img, L)
    % background subtraction
    %
    % get the min of each pixel (spatially) out of many images 
    
    % also get the mean after subtracting the min scale each image by that
    % compensates for spatial variation 
    % this is not implemented right now, I won't because it seems the 
    % current setup does not have enough images per chip for this work

    nRows = floor(size(img,1)/L);
    nCols = floor(size(img,2)/L);
    nChannels = size(img,3);
    
    for n = 1:nRows
        for m = 1:nCols
            for ci = 1:nChannels
                ymin = 1+(n-1)*L;
                ymax = min(n*L, size(img,1));
                xmin = 1+(m-1)*L;
                xmax = min(m*L, size(img,2));
                bg(:,:,ci) = min(bg(:,:,ci), img(ymin:ymax,xmin:xmax,ci));
            end
        end
    end
end