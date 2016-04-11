clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

%%
dataDir = '/Volumes/IdseData/160406_8well_sox17smad4_fixed';
channelLabel = {'Sox17','Eomes','H2B/Smad4','H2B/Sox17'};

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

for i = 18%%:numel(listing)
   
    filename = listing(i).name;
    [~, barefname] = fileparts(filename);
    metafilename = fullfile(dataDir,[barefname '_metadata.mat']);
    load(metafilename);
    
    fname = fullfile(dataDir,[barefname '_c3_Simple Segmentation.h5']);
    segGFP = separateFusedNuclei((squeeze(h5read(fname, '/exported_data')) == 2)');
    
    fname = fullfile(dataDir,[barefname '_c4_Simple Segmentation.h5']);
    segRFP = separateFusedNuclei((squeeze(h5read(fname, '/exported_data')) == 2)');
    
    img = readStack(meta.filename);
end

%%
figure, imshow(imadjust(mat2gray(img(:,:,4))))

%%
img = squeeze(img);
imgOverlay = cat(3, imadjust(mat2gray(img(:,:,4))), segGFP, segRFP);
figure, imshow(imgOverlay,[])




