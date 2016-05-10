function [MIPtot, MIPidxtot] = makeMIP_Andor(inputdir, position, channel, outputdir, saveidx)

    % turn off annoying BF warning
    warning('off','BF:lowJavaMemory');
    
    % read metadata
    meta = MetadataAndor(inputdir);
    filenameFormat = meta.filename;
    s = strsplit(meta.filename,'_');
    barefname = s{1};
    
    % some input checking
    if ~exist('saveidx','var')
        saveidx = false;
    end
    if ~exist('outputdir','var')
        outputdir = inputdir;
    end
    if ~exist(outputdir,'dir')
        mkdir(outputdir);
    end
    ci = channel;
    if ci+1 > meta.nChannels 
        error('channel exceeds number of channels');
    end
    if isempty(strfind(filenameFormat,'_w'))
        error('export files with wavelength separated');
    end
    pi = position;
    if pi+1 > meta.nPositions
        error('position exceeds number of positions');
    end
 
    MIP = {};
    MIPidx = {};

    fileTmax = ceil(meta.nTime/meta.tPerFile) - 1;
    for ti = 0:fileTmax

        % read the data 
        % (if all time points are in one file there is no _t part)
        if fileTmax > 0 
            fname = sprintf(filenameFormat, pi, ti, ci);
        else
            fname = sprintf(filenameFormat, pi, ci);
        end
        tic
        stack = readStack(fullfile(inputdir,fname));
        toc

        % make MIP
        [MIP{ti+1},MIPidx{ti+1}] = max(stack(:,:,:,:,:),[],3);
        MIPidx{ti+1} = uint16(MIPidx{ti+1});
        clear stack;
    end

    % combine time series from different files
    MIPtot = cat(5,MIP{:});
    MIPidxtot = cat(5,MIPidx{:});
    clear MIP MIPidx;

    % save result
    %-------------

    % MIP
    fname = fullfile(outputdir, sprintf([barefname '_MIP_p%.4d_w%.4d.tif'],pi,ci));
    if exist(fname,'file')
        delete(fname);
    end
    % I HAD PROBLEMS WITH BFSAVE
    %bfsave(MIPtot,fname, 'dimensionOrder', 'XYZCT');
    MIPtot = squeeze(MIPtot);
    imwrite(MIPtot(:,:,1), fname);
    for i = 2:size(MIPtot,3)
        imwrite(MIPtot(:,:,i), fname,'WriteMode','Append');
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
    
    warning('on','BF:lowJavaMemory');
end