classdef metadata
    % a class to create a uniform metadata interface

    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------

    properties

        filename
        
        xres                % x-resolution
        yres                % y-resolution
        xSize
        ySize

        nZslices

        nChannels           % number of channels
        channelNames        % cell array of channel names 
                            % when read by bioformats: detector name or dye
        excitationWavelength 

        channelLabel        % names given to channels by hand, 
                            % e.g. labeled protein

        nTime               % number of time points
        timeInterval        

        nPositions
        montageOverlap      % percent overlap of montage
        montageGridSize     % grid size n by m locations
        XYZ                 % position coordinates

        raw                 % store unprocessed metadata if necessary
    end
    
    methods
        
        function this = metadata(filename)

            if nargin == 1
                this = this.read(filename);
            end
            
            this.filename = filename;
        end
        
        function this = read(this, filename)
            % read metadata from file using bioformats
            %
            % read(filename)
            
            r = bfGetReader(filename);
            omeMeta = r.getMetadataStore();

            this.xSize = omeMeta.getPixelsSizeX(0).getValue();
            this.ySize = omeMeta.getPixelsSizeY(0).getValue();

            omeMeta = r.getMetadataStore();
            this.xres = double(omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROM));
            this.yres = double(omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROM));

            this.nChannels = r.getSizeC();
            this.channelNames = {};
            for ci = 1:this.nChannels 
                this.channelNames{ci} = char(omeMeta.getChannelName(0,ci-1));
            end
            
            this.excitationWavelength = {};
            for ci = 1:this.nChannels 
                this.excitationWavelength {ci} = round(10^3*double(omeMeta.getChannelExcitationWavelength(0,ci-1).value(ome.units.UNITS.MICROM)));
            end

            this.nZslices = r.getSizeZ();
            this.nTime = r.getSizeT();

            %omeMeta.getPixelsType(0);

            this.raw = char(omeMeta.dumpXML());
        end
        
        function save(this)
            % save this object to a mat file
            %
            % save()
            %
            % e.g. raw filename is 1.oib -> stores metadata in same place
            % under 1_metadata.mat

            meta = this;
            
            [datadir,barefname] = fileparts(this.filename);
            metafname = fullfile(datadir,[barefname '_metadata']);
            save(metafname, 'meta');
        end
    end
end