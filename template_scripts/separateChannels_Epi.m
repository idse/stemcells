%--------------------------------------------------------------------------
% script to save separate channels of btf snapshot
% for the purpose of segmentation in Ilastik
%--------------------------------------------------------------------------

% change this to the location of the code on your computer:
codeLocation = '/Users/idse/repos/Warmflash/stemcells';

% change this to the location of the data on your computer:
dataDir = '/Volumes/IdseData/160416-EndDIff-Sox17nS4cells-epi';

% change this to the channels that you want to segment (nuclei)
channels = [3 4];

%% preprocessing

warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

listing = dir(fullfile(dataDir,'*btf'));

for i = 1:numel(listing)
    
    filename = listing(i).name;
    [~, barefname] = fileparts(filename);

    for ci = channels
        img = imread(fullfile(dataDir,filename), ci);
        channelFilename = fullfile(dataDir,[barefname '_c' num2str(ci) '.tif']);
        imwrite(img, channelFilename); 
    end
end
