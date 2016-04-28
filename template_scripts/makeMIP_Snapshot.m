clear all; close all;
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
 
addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

%%

dataDir = '/Volumes/IdseData/160315_smad2';
MIPchannels = 1:2;%3:4;
saveidx = [false true]; 
previewChannels = [1 2 3 5];

%% preprocessing

MIPdir = fullfile(dataDir,'MIP');
previewDir = fullfile(dataDir,'preview');

if ~exist(MIPdir,'dir')
    mkdir(MIPdir);
end
if ~exist(previewDir,'dir')
    mkdir(previewDir);
end

listing = dir(fullfile(dataDir,'*oib'));
if isempty(listing)
    listing = dir(fullfile(dataDir,'*tif'));
end

for i = 1:numel(listing)
    
    filename = listing(i).name;
    [~, barefname] = fileparts(filename);
    
    img = readStack(fullfile(dataDir,filename));

    for cii = 1:numel(MIPchannels)
        
        ci = MIPchannels(cii);
        
        [MIP, MIPidx] = max(img(:,:,:,ci),[],3);
        
        channelFilename = fullfile(MIPdir,[barefname '_MIP_c' num2str(ci) '.tif']);
        imwrite(MIP, channelFilename); 
        
        if saveidx(cii)
            channelFilename = fullfile(MIPdir,[barefname '_MIPidx_c' num2str(ci) '.tif']);
            imwrite(MIPidx, channelFilename); 
        end
    end
    
    for cii = 1:numel(previewChannels)
        
        ci = previewChannels(cii);
        
        [MIP, MIPidx] = max(img(:,:,:,ci),[],3);
        
        channelFilename = fullfile(previewDir,[barefname '_MIP_c' num2str(ci) '.tif']);
        imwrite(imadjust(mat2gray(MIP)), channelFilename); 
    end
end
