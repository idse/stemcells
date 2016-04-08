function meta = readMeta_Andor(dataDir)
    % meta = readMeta_Andor(dataDir)
    %
    % dataDir: dir containing data and .txt metadata file
    %
    % meta: structure with fields (when applicable)
    %
    %   xres:   x-resolution
    %   yres:   y-resolution
    %   xSize:  IMPLEMENT
    %   ySize:  IMPLEMENT
    %
    %   nZslices
    %
    %   nChannels:          number of channels
    %   channelNames:       IMPLEMENT cell array of channel names in order
    %
    %   nTime:              number of time points
    %   timeInterval
    %
    %   nPositions
    %   montageOverlap:     percent overlap of montage
    %   montageGridSize:    grid size n by m locations
    %   XYZ:                position coordinates
    %
    %   Andor only file info:
    %
    %   fnameFormat         filename format as for sprintf
    %   bareFileName        filename without _suffixes
    %   tPerFile            number of time points per file
    
    %----------------------
    % Idse Heemskerk, 2016
    %----------------------

    rawMeta = struct();
    meta = struct();

    % get the image filename formate
    listing = dir(fullfile(dataDir,'*.tif'));
    filename = listing(1).name;
    [s,matches] = strsplit(filename,{'_.[0-9]{4}'},...
                           'DelimiterType','RegularExpression',...
                           'CollapseDelimiters',false);
    for i = 1:numel(matches)
        matches{i} = [matches{i}(1:2) '%.4d'];
    end
    meta.bareFileName = s{1};
    meta.fnameFormat = [s{1} matches{:} s{end}];
    
    % open the meta data file
    listing = dir(fullfile(dataDir,'*.txt'));
    filename = fullfile(dataDir, listing(1).name);
    
    fid = fopen(filename);
    
    if fid == -1
        error('file not found');
    end

    % move to the first line of the heads
    for i = 1:3
        tline = fgets(fid);
    end

    % scan through the header
    while ischar(tline) 
        k = strfind(tline,' : ');
        if isempty(k)
            break;
        else
            val = strtrim(tline(k + 3:end));
            
            if ~isnan(str2double(val))
                rawMeta.(tline(1:k-1)) = str2double(val);
            else
                rawMeta.(tline(1:k-1)) = val;
            end
        end
        tline = fgets(fid);
    end
    
    % first make a field tPerFile to keep file information together 
    if isfield(rawMeta,'Time')
        meta.tPerFile = NaN;
    end
    
    % resolution
    k = strfind(rawMeta.x,'*');
    meta.xres = str2double(rawMeta.x(k+2:k+9));
    meta.yres = str2double(rawMeta.y(k+2:k+9));
    
    % read out number of slices
    nZslices = 1;
    if isfield(rawMeta,'Z')
        nZslices = rawMeta.Z;
    end
    meta.nZslices = nZslices;

    % read out number of channels
    nChannels = 1;
    if isfield(rawMeta,'Wavelength')
        nChannels = rawMeta.Wavelength;
    end
    meta.nChannels = nChannels;
    
    % read out number of time points
    nTimePts = 1;
    if isfield(rawMeta,'Time')
        nTimePts = rawMeta.Time;
    end
    meta.nTime = nTimePts;
    
    % read out number of positions
    nPositions = 1;
    if isfield(rawMeta,'Montage')
        nPositions = rawMeta.Montage;
    elseif isfield(rawMeta,'XY')
        nPositions = rawMeta.XY;
    end
    meta.nPositions = nPositions;
    
    % read time interval
    if nTimePts > 0 
        while ischar(tline) 
            tline = fgets(fid);
            k = strfind(tline, 'Repeat T');
            if ~isempty(k)
                break;
            end
        end
        s = strsplit(tline,{'(',')'});
        meta.timeInterval = s{2};
    end
    
    % read positions
    if nPositions > 1
        
        
        % montage specific stuff
        if isfield(rawMeta,'Montage')
            
            % montage grid size
            while ischar(tline) 
                tline = fgets(fid);
                k = strfind(tline, 'Montage Positions');
                if ~isempty(k)
                    break;
                end
            end
            s = strsplit(tline,' ');
            si = find(strcmp(s,'by'));
            meta.montageGridSize = [str2double(s{si-1}) str2double(s{si+1})];
            
            % overlap for montage
            while ischar(tline) 
                tline = fgets(fid);
                k = strfind(tline, 'Overlap');
                if ~isempty(k)
                    break;
                end
            end
            overlap = str2double(strtrim(tline(k+8:end)));
            meta.montageOverlap = overlap;
        end
        
        % move on to XYZScan to read positions
        while ischar(tline) 
            if ~isempty(strfind(tline,'[Region Info (Fields)]'))
                break;
            end
            tline = fgets(fid);
        end
        for i = 1:3
            tline = fgets(fid);
        end

        % read positions
        XYZ = zeros([nPositions 3]);
        for i = 1:nPositions
            s = strsplit(tline,'\t');
            XYZ(i,1) = str2double(s{end-1});
            XYZ(i,2) = str2double(s{end});
            tline = fgets(fid);
            s = strsplit(tline,'\t');
            XYZ(i,3) = str2double(s{1});
            tline = fgets(fid);
        end
        meta.XYZ = XYZ;
    end
    
    fclose(fid);
    
    % number of time points per file
    %---------------------------------
    if isfield(rawMeta,'Time')
        
        pathstr = fileparts(filename);
        
        % determine the number of files the time series is split up into
        listing = dir(pathstr);
        info = [];
        fileTmax = 0;
        for i = 1:numel(listing)
            k = strfind(listing(i).name,'_t');
            if ~isempty(k)
                fileTmax = max(fileTmax, str2double(listing(i).name(k+2:k+5)));
                if isempty(info)
                    info = imfinfo(fullfile(pathstr,listing(i).name));
                end
            end
        end
        % if files don't have a _t label nothing was found
        if isempty(info)
            meta.tPerFile = meta.nTime;
        else
            % time points per file
            meta.tPerFile = numel(info)/meta.nZslices;
        end
    end
end