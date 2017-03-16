function batchMIP_lsm(inputdir, outputdir, channels, saveidx, tmax, positions)

    if ~exist('positions','var')
        tracks = dir(fullfile(inputdir,'Track*')); 
        positions = 1:numel(tracks);
    end

    if ~exist('tmax','var')
        tracks = dir(fullfile(inputdir,'Track*'));
        files = dir(fullfile(inputdir, tracks(1).name, '*oif'));
        tmax = numel(files);
    end

    for ci = 1:numel(channels)
        for pi = positions
            disp(['processing position ' num2str(pi) ', channel ' num2str(channels(ci))]);
            makeMIP_lsm(inputdir, pi-1, channels(ci), outputdir, saveidx(ci), tmax);
        end
    end
end