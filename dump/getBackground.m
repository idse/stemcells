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

% BACKGROUND FROM INTENSITY DISTRIBUTION
% %%
% minI = min(S4Stitched(S4Stitched>0));
% maxI = 5000;
% figure,imshow(S4Stitched - 811,[0 maxI-811])
% 
% %%
% figure,
% bins = 1:5:maxI;
% n = histc(double(S4Stitched(:)-minI),bins);
% bar(bins,n);
% axis([bins(1) bins(end) 0 10^5]);
% 
% %%
% %[~,i] = max(conv(n,GaussD(5,1)','same')')
% bla = conv(n,GaussD(5,2)','same')';
% [~,i] = min(bla);
% plot(bla)
% bins(i)
% 
% %%
% x = medfilt2(S4Stitched,[3 3]);
% n = histc(double(x(:)-minI),bins);
% bar(bins,n);
% axis([bins(1) bins(end) 0 10^5]);
% bla = conv(n,GaussD(5,2)','same')';
% plot(bla)
% 
% %%
% imshow(x-816,[0 5000])