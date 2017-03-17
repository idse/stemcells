function [MIP, MIPidxtot] = makeMIP_lsm(inputdir, position, channel, outputdir, saveidx, tmax)

    % positions are called Track and stored in separate directories
    tracks = dir(fullfile(inputdir,'Track*'));
    files = dir(fullfile(inputdir, tracks(1).name, '*oif'));
    fname = fullfile(inputdir, tracks(1).name, files(1).name);
    
    meta = Metadata(fname);
    meta.nTime = numel(files);

    % some input checking
    if ~exist('saveidx','var')
        saveidx = false;
    end
    if ~exist('outputdir','var') || isempty(outputdir)
        outputdir = fullfile(dataDir,'MIP');
    end
    if ~exist(outputdir,'dir')
        mkdir(outputdir);
    end
    ci = channel;
    if ~exist('tmax','var')
        tmax = meta.nTime;
    end

    ci = channel;
    pi = position;

    files = dir(fullfile(inputdir, tracks(pi+1).name, '*oif'));
    s = strsplit(files(1).name,'_');
    barefname = s{1};
    
    MIP = zeros([meta.ySize meta.xSize tmax], 'uint16');
    MIPidxtot = zeros([meta.ySize meta.xSize tmax], 'uint16');

    for ti = 1:tmax

        fname = fullfile(inputdir, tracks(pi+1).name, sprintf([barefname '_%.2d.oif'], ti));
        r = bfGetReader(fname);

        if r.getSizeZ() > 1
            im = zeros([r.getSizeY() r.getSizeX() r.getSizeZ()]);
            for zi = 1:r.getSizeZ()
                im(:,:,zi) = bfGetPlane(r, r.getIndex(zi-1,ci,0)+1);
            end
            [MIP(:,:,ti), MIPidxtot(:,:,ti)] = max(im,[],3);
        else
            MIP(:,:,ti) = bfGetPlane(r, r.getIndex(0,ci,0)+1);
        end
        
        fprintf('.');
        if mod(ti,60)==0
            fprintf('\n');
        end
    end
    fprintf('\n');

    % save result
    %-------------

    % MIP
    s = strsplit(inputdir,filesep);
    barefname = s{end};
    fname = fullfile(outputdir, sprintf([barefname '_MIP_p%.4d_w%.4d.tif'],pi,ci));
    if exist(fname,'file')
        delete(fname);
    end

    MIP = squeeze(MIP);
    imwrite(MIP(:,:,1), fname);
    for i = 2:size(MIP,3)
        imwrite(MIP(:,:,i), fname,'WriteMode','Append');
    end

    % MIP index
    if saveidx
        fname = fullfile(outputdir, sprintf([barefname '_MIPidx_p%.4d_w%.4d.tif'],pi,ci));
        if exist(fname,'file')
            delete(fname);
        end
        %bfsave(MIPidxtot,fname, 'dimensionOrder', 'XYZCT');
        MIPidxtot = squeeze(MIPidxtot);
        imwrite(MIPidxtot(:,:,1), fname);
        for i = 2:size(MIPidxtot,3)
            imwrite(MIPidxtot(:,:,i), fname, 'WriteMode','Append');
        end
    end
end