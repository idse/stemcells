% produces a binary image showing areas that are separated from the rest of the
% image by a relatively sharp change in grayscale value.
function seg = imseg(im,diskSize,stdDev,minArea)

% find difference between min and max value around each pixel
diff = imdilate(im,strel('disk',diskSize)) - imerode(im,strel('disk',diskSize));

% find where the diff is unusually large
segx = diff > (mean(mean(diff)) + stdDev * std(diff(:)));

% fill in the bright spots by different padding combinations
segx_a = padarray(segx,[1 1],1,'pre');
segx_a_filled = imfill(segx_a,'holes');
segx_a_filled = segx_a_filled(2:end,2:end);

segx_b = padarray(padarray(segx,[1 0],1,'pre'),[0 1],1,'post');
segx_b_filled = imfill(segx_b,'holes');
segx_b_filled = segx_b_filled(2:end,1:end-1);

segx_c = padarray(segx,[1 1],1,'post');
segx_c_filled = imfill(segx_c,'holes');
segx_c_filled = segx_c_filled(1:end-1,1:end-1);

segx_d = padarray(padarray(segx,[1 0],1,'post'),[0 1],1,'pre');
segx_d_filled = imfill(segx_d,'holes');
segx_d_filled = segx_d_filled(1:end-1,2:end);

% put them all together
segx_filled = segx_a_filled | segx_b_filled | segx_c_filled | segx_d_filled;
seg = segx_filled;
%seg = diff > (mean(mean(diff)) + stdDev * std(diff(:))); % no filling

% take out little spots
seg = bwareaopen(seg,minArea);

% overcompensate a little
seg = imdilate(seg,strel('disk',1));

return