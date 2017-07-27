classdef Metadata
    % a class to create a uniform metadata interface

    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------

    properties

        filename
        
        %tPerFile
        
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
        
        % to keep track of multi-well experiments
        nWells              % number of well
        posPerCondition     % positions per well
        conditions          % cell array of condition labels
    end
    
    methods
        
        function this = Metadata(filename)

            if nargin == 1
                this = this.read(filename);
                this.filename = filename;
            end
        end
        
        function this = read(this, filename)
            % read metadata from file using bioformats
            %
            % read(filename)
            
            r = bfGetReader(filename);
            omeMeta = r.getMetadataStore();

            this.xSize = omeMeta.getPixelsSizeX(0).getValue();
            this.ySize = omeMeta.getPixelsSizeY(0).getValue();

            this.xres = double(omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROM));
            this.yres = double(omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROM));

            dt = omeMeta.getPixelsTimeIncrement(0);
            if ~isempty(dt)
                this.timeInterval = double(dt.value);
            else
                this.timeInterval = [];
            end
            
            this.nChannels = r.getSizeC();
            this.channelNames = {};
            for ci = 1:this.nChannels 
                this.channelNames{ci} = char(omeMeta.getChannelName(0,ci-1));
            end
            
            this.excitationWavelength = {};
            for ci = 1:this.nChannels 
                lambda = omeMeta.getChannelExcitationWavelength(0,ci-1);
                if ~isempty(lambda)
                    this.excitationWavelength{ci} = round(10^3*double(lambda.value(ome.units.UNITS.MICROM)));
                end
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

            [datadir,barefname] = fileparts(this.filename);
            metafname = fullfile(datadir,[barefname '_metadata']);
            save(metafname, 'this');
        end
        
        function displayPositions(this)
            % display positions 
            %
            % displayPositions();
            %
            % positions are center of field of view?
            
            XYZ = this.XYZ;

            %scatter(XYZ(:,1), XYZ(:,2))

            for i = 1:size(XYZ,1)
                text(XYZ(i,1), XYZ(i,2),num2str(i))
                w = 1024*this.xres;
                h = 1024*this.yres;
                rectangle('Position',[XYZ(i,1)-w/2,XYZ(i,2)-h/2,w,h])
            end
            axis([min(XYZ(:,1))-w max(XYZ(:,1))+w min(XYZ(:,2))-h max(XYZ(:,2))+h])
            axis equal
            axis off

            % XYZmean = mean(XYZ);
            % hold on
            % scatter(XYZmean(:,1), XYZmean(:,2),'r')
            % hold off
        end
    end
end