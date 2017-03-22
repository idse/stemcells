function img = readStack(fullfname, channels) 
    % read the data from a single multichannel stack
            
    % load the Bio-Formats library into the MATLAB environment
    autoloadBioFormats = 1;
    status = bfCheckJavaPath(autoloadBioFormats);
    assert(status, ['Missing Bio-Formats library. Either add loci_tools.jar '...
        'to the static Java path or add it to the Matlab path.']);

    % initialize logging
    loci.common.DebugTools.enableLogging('INFO');
    
    % create bioformats reader for file
    disp(fullfname);
    
%     r = bfGetReader(fullfname);
% 
%     imclass = class(bfGetPlane(r, 1));
%     
%     stackSize = [r.getSizeY(), r.getSizeX(), r.getSizeZ(), r.getSizeC(), r.getSizeT()];
%     data = zeros(stackSize, imclass);
%     
%     for i = 1:r.getImageCount()
% 
%         ZCTidx = r.getZCTCoords(i-1) + 1;
%         
%         fprintf('.');
%         if rem(i,80) == 0
%             fprintf('\n');
%         end
% 
%         data(:,:, ZCTidx(1), ZCTidx(2), ZCTidx(3)) = bfGetPlane(r, i);
%     end
%     fprintf('\n');
% 
%     nChannels = r.getSizeC();
%     r.close();

    [~,~,ext] = fileparts(fullfname);

    if strcmp(ext,'.tif') || strcmp(ext,'.btf')

        info = imfinfo(fname);
        w = info.Width;
        h = info.Height;

        img = zeros([h w numel(channels)],'uint16');
        for cii = 1:numel(channels)
            img(:,:,cii) = imread(fullfname, channels(cii));
        end

        if exist('time','var') && time > 1
            error('todo : include reading for dynamic not Andor or epi');
        end

    elseif strcmp(ext,'.vsi') || strcmp(ext, '.oif')

        r = bfGetReader(fullfname);
        img = zeros([r.getSizeY() r.getSizeX() numel(channels) r.getSizeZ()], 'uint16');
        for cii = 1:numel(channels)
            for zi = 1:r.getSizeZ()
                time = 1; % intended for snapshots
                img(:,:,cii,zi) = bfGetPlane(r, r.getIndex(zi-1,channels(cii)-1,time-1)+1);
            end
        end
        r.close();
        img = squeeze(img);
    end
end
