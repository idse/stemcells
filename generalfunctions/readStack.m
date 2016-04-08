function [data nChannels] = readStack(fullfname) 
    % read the data from a single multichannel stack

    % a way to do it without bioformats:
    %         tic
    %         disp(['reading: ' fname]);
    %         warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning');
    %         stack = zeros([1024 1024 meta.nZslices meta.tPerFile]);
    %         imax = meta.nZslices*meta.tPerFile;
    %         for i = 0:imax-1
    %             zi = rem(i,meta.nZslices)+1;
    %             tj = ceil((i+1)/meta.nZslices);
    %             stack(:,:,zi,tj) = imread(fullfile(inputdir,fname),i+1);
    %         end
    %         warning('on', 'MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning');
    %         toc
            
    % load the Bio-Formats library into the MATLAB environment
    autoloadBioFormats = 1;
    status = bfCheckJavaPath(autoloadBioFormats);
    assert(status, ['Missing Bio-Formats library. Either add loci_tools.jar '...
        'to the static Java path or add it to the Matlab path.']);

    % initialize logging
    loci.common.DebugTools.enableLogging('INFO');
    
    % create bioformats reader for file
    disp(fullfname);
    
    r = bfGetReader(fullfname);

    imclass = class(bfGetPlane(r, 1));
    
    stackSize = [r.getSizeY(), r.getSizeX(), r.getSizeZ(), r.getSizeC(), r.getSizeT()];
    data = zeros(stackSize, imclass);
    
    for i = 1:r.getImageCount()

        ZCTidx = r.getZCTCoords(i-1) + 1;
        
        fprintf('.');
        if rem(i,80) == 0
            fprintf('\n');
        end

        data(:,:, ZCTidx(1), ZCTidx(2), ZCTidx(3)) = bfGetPlane(r, i);
    end
    fprintf('\n');

    nChannels = r.getSizeC();
    r.close();
end
