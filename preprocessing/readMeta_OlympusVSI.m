function meta = readMeta_OlympusVSI(vsifile)
    % meta = readMeta_Andor(dataDir)
    %
    % vsifile:  full path to Olympus .vsi file
    %
    % meta: structure with fields (when applicable)
    %
    %   xres:   x-resolution
    %   yres:   y-resolution
    %   xSize
    %   ySize
    %
    %   nZslices
    %
    %   nChannels:          number of channels
    %   channelNames:       cell array of channel names in order
    %
    %   nTime:              number of time points
    %   timeInterval        IMPLEMENT
    %
    %   nPositions
    %   montageOverlap:     percent overlap of montage
    %   montageGridSize:    grid size n by m locations
    %   XYZ:                position coordinates
    %
    
    %----------------------
    % Idse Heemskerk, 2016
    %----------------------
    
    meta = struct();
    
    r = bfGetReader(vsifile);
    omeMeta = r.getMetadataStore();
    
    meta.xSize = omeMeta.getPixelsSizeX(0).getValue();
    meta.ySize = omeMeta.getPixelsSizeY(0).getValue();

    omeMeta = r.getMetadataStore();
    meta.xres = double(omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROM));
    meta.yres = double(omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROM));

    meta.nChannels = r.getSizeC();
    meta.channelNames = {};
    for ci = 1:meta.nChannels 
        meta.channelNames{ci} = char(omeMeta.getChannelName(0,ci-1));
    end
    
	meta.nZslices = r.getSizeZ();
    meta.nTime = r.getSizeT();
    
    %omeMeta.getPixelsType(0);
    
    meta.omeXML = char(omeMeta.dumpXML());
end