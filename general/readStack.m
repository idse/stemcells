function [data meta] = readStack(fullfname) 
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

    meta = struct();
    meta.nChannels = r.getSizeC();
    
    omeMeta = r.getMetadataStore();
    meta.xres = double(omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROM));
    meta.yres = double(omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROM));

    r.close();
end
