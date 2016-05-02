clear all; close all;
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
 
addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

%%

%dataDir = '/Volumes/IdseData/160315_smad2';
dataDir = '/Volumes/IdseData/160416-Pagl vs RPMI epi';
MIPchannels = [1 2 3 4];%1:2;%3:4;
saveChunk = true;
saveidx = [false false false false]; 

previewChannels = [];

%% preprocessing

MIPdir = fullfile(dataDir,'MIP');
previewDir = fullfile(dataDir,'preview');


if ~exist(MIPdir,'dir') && ~isempty(MIPchannels)
    mkdir(MIPdir);
end
if ~exist(previewDir,'dir')
    mkdir(previewDir);
end


listing = dir(fullfile(dataDir,'*oib'));
if isempty(listing)
    listing = dir(fullfile(dataDir,'*tif'));
end
if isempty(listing)
    listing = dir(fullfile(dataDir,'*btf'));
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
        
        if saveChunk
            chunk = MIP(2048:2*2048, 2048:2*2048);
            channelFilename = fullfile(MIPdir,[barefname '_chunk_c' num2str(ci) '.tif']);
            imwrite(chunk, channelFilename); 
        end
    end

    for cii = 1:numel(previewChannels)
        
        ci = previewChannels(cii);
        
        [MIP, MIPidx] = max(img(:,:,:,ci),[],3);
        
        channelFilename = fullfile(previewDir,[barefname '_c' num2str(ci) '.tif']);
        imwrite(imadjust(mat2gray(MIP)), channelFilename); 
    end
end
