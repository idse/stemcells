classdef DynamicPositionAndor < Position
    % Data class to store cell data in a field of view

    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------
    
    methods

        function this = DynamicPosition(meta, ID)
            % constructor 
            %
            % DynamicPositionAndor(meta, ID)
            %
            % meta:     MetadataAndor object
            % ID:       position index
            
            % matlab sucks
            if nargin == 0
                return
            end
            
            this.nChannels = meta.nChannels;
            this.cellData = struct();
            
            this.nTime = meta.nTime;
            this.ID = ID;
            
            % this is a clunky way to only put the actual position in the
            % filename format but leave the rest of the %.4d pieces
            filenameFormat = meta.filename;
            barefname = sprintf(filenameFormat,this.ID-1);
            barefname = barefname(1:end-2);
            filenameFormat(1:numel(barefname)) = barefname;
            this.filename = filenameFormat;
        end

        % saving and loading
        %---------------------------------

        function img = loadImage(this, dataDir, channel, time)
            % load image of colony
            %
            % img = loadImage(dataDir, channel, time)
            %
            % dataDir:  main data directory 
            % channel:  desired channel to be loaded
            %
            % img:      loaded image
            %
            % for now assume Andor format input, can be expanded later

            % fti : time index of file, e.g. if tPerFile = 2
            % Andor times start at 0, our time index and subti at 1
            % subti : time index within file
            %
            % so time=1: fti = 0, subti = 1
            % time=30: fti = 16, subti = 1
            
            warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning');

            fname = fullfile(dataDir, sprintf(this.filename,channel-1));
            nZslices = numel(info)/this.nTime;
            ioffset = nZslices*(time-1);
            
            img = zeros([h w nZslices],'uint16');
            for i = 1:nZslices
                img(:,:,i) = imread(fname, ioffset + i);
            end
            
            %disp(['loaded image ' fname]);
        end
    end
end