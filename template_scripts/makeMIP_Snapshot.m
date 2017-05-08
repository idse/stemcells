clear all; close all;
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
 
addpath(genpath('~/Documents/Stemcells')); 

dataDir = '/Volumes/IdseData3/170227_3daysManystains/d1_10pm';

MIPchannels = 4;
saveChunk = true;
saveidx = false; 
previewChannels = 4;

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
    listing = dir(fullfile(dataDir,'*vsi'));
    meta = Metadata(fullfile(dataDir,listing(1).name));
end
if isempty(listing)
    listing = dir(fullfile(dataDir,'*tif'));
end

for i = 1:numel(listing)
    
    filename = listing(i).name;
    [~, barefname] = fileparts(filename);
    s = strsplit(barefname,'_p');
    barefname = s{1};
    
    if numel(s)>1
        posnr = s{2};
        postfix = ['_p' posnr];
    else
        postfix = [];
    end

    allChannels = unique([MIPchannels previewChannels]);
    
    for cii = 1:numel(allChannels)

        ci = allChannels(cii);
        
        disp(['channel: ' num2str(ci)]);
        img = readStack2(fullfile(dataDir,filename), ci);

        disp('making MIP');
        
        [MIP, MIPidx] = max(img,[],3);
        MIPidx = uint8(MIPidx);
        ymin = 1;
        ymax = min(ymin + 2*2048, size(MIP,1));
        xmin = 1;
        xmax = min(xmin + 2*2048, size(MIP,2));
        chunk = MIP(ymin:ymax, xmin:xmax);

        disp('saving MIP');
        
        if any(MIPchannels == ci)
           
            channelFilename = fullfile(MIPdir, sprintf([barefname '_MIP' postfix '_w%.4d.tif'], ci-1));
            imwrite(MIP, channelFilename); 

            if saveidx(MIPchannels == ci)
                channelFilename = fullfile(MIPdir, sprintf([barefname '_MIPidx' postfix '_w%.4d.tif'], ci-1));
                imwrite(MIPidx, channelFilename); 
            end
        end
        
        disp('saving preview');
        if any(previewChannels == ci)
            
            channelFilename = fullfile(previewDir, sprintf([barefname postfix '_w%.4d.tif'], ci-1));
            imwrite(imadjust(mat2gray(chunk)), channelFilename); 
        end
        
        disp('saving MIP chunck');
        if saveChunk(MIPchannels == ci)
            
            channelFilename = fullfile(MIPdir, sprintf([barefname '_chunk' postfix '_w%.4d.tif'], ci-1));
            imwrite(chunk, channelFilename); 
        end
    end
end
