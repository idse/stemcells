function batchMIP_epi(inputdir, outputdir, channels, saveidx, tmax, positions)

    % determine process numbers of vsi files
    vsifiles = dir(fullfile(inputdir,'*.vsi'));
    processNumbers = zeros([numel(vsifiles) 1],'uint16');
    for i = 1:numel(vsifiles)
        s = strsplit(vsifiles(i).name,{'_','.'});
        processNumbers(i) = uint16(str2double(s{2}));
    end
    %processnr = min(processNumbers):max(processNumbers);
    processnr = processNumbers;
    
    % read metadata
    vsifile = fullfile(inputdir,['Process_' num2str(processnr(1)) '.vsi']);
    meta = Metadata(vsifile);

    if ~exist('positions','var')
        positions = 1:numel(processnr);
    end

    if ~exist('tmax','var')
        tmax = meta.nTime;
    end

    for ci = 1:numel(channels)
        for pi = positions
            disp(['processing position ' num2str(pi) ', channel ' num2str(channels(ci))]);
            vsifile = fullfile(inputdir,['Process_' num2str(processnr(pi)) '.vsi']);
            makeMIP_epi(vsifile, pi-1, channels(ci), outputdir, saveidx(ci), tmax);
        end
    end
end