clear all; close all;
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

%%
dataDir = '/Volumes/IdseData/160406_8well_sox17smad4_fixed';
channelLabel = {'Sox17','Eomes','H2B/Smad4','H2B/Sox17'};

nChannels = numel(channelLabel);

%% preprocessing

listing = dir(fullfile(dataDir,'*oib'));

for i = 1:numel(listing)
    
    filename = listing(i).name;
    [~, barefname] = fileparts(filename);
    metafilename = fullfile(dataDir,[barefname '_metadata.mat']);
    
    if ~exist(metafilename,'file')
    
        
        filename = fullfile(dataDir,filename);

        meta = metadata(filename);
        meta.channelLabel = channelLabel;

        meta.save();

        img = readStack(meta.filename);
        
        for ci = 3:4
            channelFilename = fullfile(dataDir,[barefname '_c' num2str(ci) '.tif']);
            imwrite(max(img(:,:,:,ci),[],3), channelFilename); 
        end
    end
end

%%

listing = dir(fullfile(dataDir,'*oib'));
filename = listing(i).name;
[~, barefname] = fileparts(filename);

test = position(nChannels, filename);

img = test.loadImage(dataDir);
seg = test.loadSegmentation(dataDir);

test.extractData(dataDir);

%%
for i = 18%%:numel(listing)
   
    filename = listing(i).name;
    [~, barefname] = fileparts(filename);
    metafilename = fullfile(dataDir,[barefname '_metadata.mat']);
    load(metafilename);
end

%%
figure, imshow(imadjust(mat2gray(img(:,:,4))))

%%
img = squeeze(img);
imgOverlay = cat(3, imadjust(mat2gray(img(:,:,4))), seg(:,:,1), seg(:,:,2));
figure, imshow(imgOverlay,[])




