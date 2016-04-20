clear all; close all;
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
 
addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

%%

dataDir = '/Volumes/IdseData/160406_8well_sox17smad4_fixed';
MIPchannels = 3:4;

%% preprocessing

listing = dir(fullfile(dataDir,'*oib'));

for i = 1:numel(listing)
    
    filename = listing(i).name;
    [~, barefname] = fileparts(filename);
    
    img = readStack(meta.filename);

    for ci = MIPchannels
        channelFilename = fullfile(dataDir,[barefname '_c' num2str(ci) '.tif']);
        imwrite(max(img(:,:,:,ci),[],3), channelFilename); 
    end
end
