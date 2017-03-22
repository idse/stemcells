clear all; close all;
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
 
addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

dataDir = '/Users/Idse/data_tmp/cellfate/pulses/170222_SyringePumpAendo2';

MIPchannels = 2;
saveChunk = true;
saveidx = false; 
previewChannels = 1:2;

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
        img = readStack(fullfile(dataDir,filename), ci);

        disp('making MIP');
        
        [MIP, MIPidx] = max(img,[],3);
        MIPidx = uint8(MIPidx);
        chunk = MIP(2048:2*2048, 2048:2*2048);

        if any(MIPchannels == ci)
           
            channelFilename = fullfile(MIPdir, sprintf([barefname '_MIP' postfix '_w%.4d.tif'], ci-1));
            imwrite(MIP, channelFilename); 

            if saveidx(MIPchannels == ci)
                channelFilename = fullfile(MIPdir, sprintf([barefname '_MIPidx' postfix '_w%.4d.tif'], ci-1));
                imwrite(MIPidx, channelFilename); 
            end
        end
        
        if any(previewChannels == ci)
            
            channelFilename = fullfile(previewDir, sprintf([barefname postfix '_w%.4d.tif'], ci-1));
            imwrite(imadjust(mat2gray(chunk)), channelFilename); 
        end
        
        if saveChunk(MIPchannels == ci)
            
            channelFilename = fullfile(MIPdir, sprintf([barefname '_chunk' postfix '_w%.4d.tif'], ci-1));
            imwrite(chunk, channelFilename); 
        end
    end
end
