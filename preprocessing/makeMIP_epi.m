function [MIP, MIPidxtot] = makeMIP_epi(vsifile, position, channel, outputdir, saveidx, tmax)

    r = bfGetReader(vsifile);
    dataDir = fileparts(vsifile);

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
    if ci+1 > r.getSizeC()
        error('channel exceeds number of channels');
    end
    if ~exist('tmax','var')
        tmax = r.getSizeT();
    end

    MIPidxtot = [];
    if saveidx 
        error('TODO: implement saveidx for epi');
    end

    MIP = zeros([r.getSizeY() r.getSizeX() tmax], 'uint16');
    MIPidx = zeros([r.getSizeY() r.getSizeX() tmax], 'uint16');
    for ti = 1:tmax
        if r.getSizeZ() > 1
            im = zeros([r.getSizeY() r.getSizeX() r.getSizeZ()]);
            for zi = 1:r.getSizeZ()
                im(:,:,zi) = bfGetPlane(r, r.getIndex(zi-1,ci,ti-1)+1);
            end
            [MIP(:,:,ti), MIPidx(:,:,ti)] = max(im,[],3);
        else
            MIP(:,:,ti) = bfGetPlane(r, r.getIndex(0,ci,ti-1)+1);
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
    s = strsplit(dataDir,filesep);
    barefname = s{end};
    pi = position;
    fname = fullfile(outputdir, sprintf([barefname '_MIP_p%.4d_w%.4d.tif'],pi,ci));
    if exist(fname,'file')
        delete(fname);
    end
    
    MIP = squeeze(MIP);
    imwrite(MIP(:,:,1), fname);
    for i = 2:size(MIP,3)
        imwrite(MIP(:,:,i), fname,'WriteMode','Append');
    end
end