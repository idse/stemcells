clear all; close all;
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
 
addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

%%

dataDir = '/Volumes/IdseData/160416-EndDIff-Sox17nS4cells-epi';
MIPchannels = 3:4;

%% preprocessing

listing = dir(fullfile(dataDir,'*btf'));

for i = 1:numel(listing)
    
    filename = listing(i).name;
    [~, barefname] = fileparts(filename);

    for ci = MIPchannels
        img = imread(fullfile(dataDir,filename), ci);
        channelFilename = fullfile(dataDir,[barefname '_c' num2str(ci) '.tif']);
        imwrite(img, channelFilename); 
    end
end
