% produces a binary image showing where grayscale values are unusually high
function seg = imseg2(im,diskSize,stdDev,minArea)

% find statistics
exp = mean(mean(im));
dev = std(im(:));

% find where values are unusually large
seg = im > (exp + stdDev * dev);

% take out little spots
seg = bwareaopen(seg,minArea);

% overcompensate a little
seg = imdilate(seg,strel('disk',diskSize));

return