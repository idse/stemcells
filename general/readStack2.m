function [img, omeMeta] = readStack2(fullfname, channels, tmax) 
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

  %  [~,~,ext] = fileparts(fullfname);

%     r = bfGetReader(fullfname);
%     img = zeros([r.getSizeY() r.getSizeX() numel(channels) r.getSizeZ()], 'uint16');
%     for cii = 1:numel(channels)
%         for zi = 1:r.getSizeZ()
%             time = 1; % intended for snapshots
%             img(:,:,cii,zi) = bfGetPlane(r, r.getIndex(zi-1,channels(cii)-1,time-1)+1);
%         end
%     end
%     r.close();
%     img = squeeze(img);
    
    r = bfGetReader(fullfname);
    if ~exist('channels','var') || isempty(channels)
        channels = 1:r.getSizeC();
    end
    if ~exist('tmax','var')
        tmax = r.getSizeT();
    end
    img = zeros([r.getSizeY() r.getSizeX() numel(channels) r.getSizeZ() tmax], 'uint16');
    for ti = 1:tmax
        for cii = 1:numel(channels)
            for zi = 1:r.getSizeZ()
                img(:,:,cii,zi,ti) = bfGetPlane(r, r.getIndex(zi-1,channels(cii)-1,ti-1)+1);
            end
        end
    end
    omeMeta = r.getMetadataStore();
    r.close();
end
