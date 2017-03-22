function dataCombined = combineQPCR(dataDir, filenames, tolerance)
    % combineQPCR(dataDir, filenames, tolerance)
    
    % SOP : StepOne Plus
    data = {};
    targetsComb = {};
    samplesComb = {};
    
    for i = 1:numel(filenames)
        
        data{i} = readSOPdata3(fullfile(dataDir,filenames{i}), tolerance);
        
        targetsComb = unique(cat(1,targetsComb,data{i}.targets),'stable');
        samplesComb = unique(cat(1,samplesComb,data{i}.samples),'stable');
        
        data{i}
    end

    %data{1}.targets{3} = 'NANOG';

    % collect CT values etc
    CTmean = [];
    CTstd = [];
    Ntargets = 0;
    Nsamples = 0;
    targets = [];
    samples = [];

    if numel(targetsComb) >= numel(data{1}.targets) &&...
            numel(samplesComb) == numel(data{1}.samples)
        
        if numel(filenames) > 1
            disp('combining multiple plates with different targets, same samples');
        end
        for i = 1:numel(data)

            D = data{i};

            CTmean = cat(2, CTmean, D.CT);
            CTstd = cat(2, CTstd, D.CTstd);

            Ntargets = Ntargets + D.Ntargets;
            targets = cat(1, targets, D.targets);

            if i == 1
                Nsamples = D.Nsamples;
                samples = D.samples;
            elseif Nsamples ~= D.Nsamples
                error('samples dont match');
            end
        end

    elseif numel(targetsComb) == numel(data{1}.targets) &&...
            numel(samplesComb) > numel(data{1}.samples)

        disp('combining multiple plates with same targets, different samples');

        for i = 1:numel(data)

            D = data{i};

            CTmean = cat(1, CTmean, D.CT);
            CTstd = cat(1, CTstd, D.CTstd);

            Nsamples = Nsamples + D.Nsamples;
            samples = cat(1, samples, D.samples);

            if i == 1
                Ntargets = D.Ntargets;
                targets = D.targets;
            elseif Ntargets ~= D.Ntargets
                error('samples dont match');
            end
        end
    end

    dataCombined = struct(...
        'targets',{targets},'samples',{samples},'CTmean',CTmean,'CTstd',CTstd,...
        'Nsamples', Nsamples, 'Ntargets', Ntargets);
end