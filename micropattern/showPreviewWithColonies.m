function showPreviewWithColonies(dataDir,imagetype)
% Shows image with overlay of identified colonies. 
% dataDir is a directory produced by the processVsi function
% imagetype must be DAPI or RGB
%
% see also: processVsi


if ~exist('imagetype','var')
    imagetype = 'DAPI';
end

if strcmpi(imagetype,'DAPI')
    preview = imread(fullfile(dataDir,'previewDAPI.tif'));
    preview = imadjust(preview);
elseif strcmpi(imagetype,'RGB')
    preview = imread(fullfile(dataDir,'previewRGB.tif'));
    for ii = 1:size(preview,3)
        preview(:,:,ii) = imadjust(preview(:,:,ii));
    end
else
    disp('Image type must be DAPI or RGB');
end

load(fullfile(dataDir,'colonies.mat'));
load(fullfile(dataDir,'metaData.mat'));
figure,
imshow(preview);
hold on
CM = cat(1,colonies.center);
scale = size(preview,1)/meta.ySize;
CM(:,2) = CM(:,2)*scale;
CM(:,1) = CM(:,1)*scale;
radius = cat(1,colonies.radiusPixel)*scale;    
viscircles(CM,radius,'LineWidth',1)
for i = 1:size(CM,1)
   text(CM(i,1),CM(i,2), num2str(i),'Color','red','BackgroundColor','white',...
       'Margin',1,'FontSize',8,'HorizontalAlignment','center'); 
end