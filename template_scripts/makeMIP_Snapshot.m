clear all; close all;
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
 
addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

dataDir = '/Users/idse/data_tmp/170311_Smad2Dose';
%dataDir = '/Users/idse/data_tmp/160901_smad2';
%dataDir = '/Volumes/IdseData3/161009_Smad2FISH';
%dataDir = '/Users/idse/data_tmp/161009_Smad2FISH';

MIPchannels = 1:3;%1:2;%3:4;
saveChunk = false;
%saveidx = [true false]; 
saveidx = [false true false]; 

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
    s = strsplit(barefname,'_p');
    barefname = s{1};
    posnr = s{2};
    
    img = readStack(fullfile(dataDir,filename));

    disp('making MIPs');
    for cii = 1:numel(MIPchannels)
        
        ci = MIPchannels(cii);
        disp(['channel: ' num2str(ci)]);
        
        [MIP, MIPidx] = max(img(:,:,:,ci),[],3);
        MIPidx = uint8(MIPidx);
        
        channelFilename = fullfile(MIPdir, sprintf([barefname '_MIP_p' posnr '_w%.4d.tif'], ci-1));
        imwrite(MIP, channelFilename); 
        
        if saveidx(cii)
            channelFilename = fullfile(MIPdir, sprintf([barefname '_MIPidx_p' posnr '_w%.4d.tif'], ci-1));
            imwrite(MIPidx, channelFilename); 
        end
        
        if saveChunk
            chunk = MIP(2048:2*2048, 2048:2*2048);
            channelFilename = fullfile(MIPdir, sprintf([barefname '_chunk_p' posnr '_w%.4d.tif'], ci-1));
            imwrite(chunk, channelFilename); 
        end
    end

    disp('making previews');
    for cii = 1:numel(previewChannels)
        
        ci = previewChannels(cii);
        
        [MIP, MIPidx] = max(img(:,:,:,ci),[],3);
        
        channelFilename = fullfile(previewDir, sprintf([barefname '_p' posnr '_w%.4d.tif'], ci-1));
        imwrite(imadjust(mat2gray(MIP)), channelFilename); 
    end
end
