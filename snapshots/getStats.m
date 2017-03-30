function stats = getStats(dataDir, filelist, normalizeChannel)

    if exist('normalizeChannel','var')
        normalize = true;
    else
        normalize = false;
    end
    
    allData = {};
    nucLevel = {};
    nucLevelAll = [];

    for fi = 1:numel(filelist)

        vsifile = fullfile(dataDir, filelist{fi});
        [~,barefname,~] = fileparts(vsifile);

        load(fullfile(dataDir,'results', [barefname '_positions.mat']));
        allData{fi} = P;

        if normalize
            N = P.cellData.nucLevel(:,normalizeChannel);
            nucLevel{fi} = bsxfun(@rdivide, P.cellData.nucLevel, N);
        else
            nucLevel{fi} = P.cellData.nucLevel;
        end
        nucLevelAll = cat(1, nucLevelAll, nucLevel{fi});
    end

    % determine limits for displaying intensities
    tol = 0.001;
    lim = {};
    for i = 1:P.nChannels
        A = nucLevelAll(:,i);
        lim{i} = stretchlim(mat2gray(A),tol)*(max(A) - min(A)) + min(A);
    end

    % histogram of different channels
    [bins, hist] = makeHistograms(nucLevel, lim);

    stats = struct('histograms',{hist},'bins',{bins},...
        'lim',{lim}, 'nucLevel', {nucLevel});
end