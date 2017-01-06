function data = readSOPdata3(filename)

    fid = fopen(filename,'r','n','UTF-8');

    tline = fgetl(fid);

    firstlinefound = 0;
    q = 1;
    while(ischar(tline))

        if firstlinefound == 0
             if strfind(tline,'Sample')
                 firstlinefound = 1;
                 rawdata = strsplit(tline,'\t');
             end
        else
            if ~isempty(tline)
                rawdata = cat(1, rawdata, strsplit(tline,'\t','CollapseDelimiters',false));
                q = q + 1;
            else
                break;
            end
        end
        tline = fgetl(fid);
        tline = strrep(tline,',',', ');
    end
   
    Scoli = find(~cellfun(@isempty,strfind(rawdata(1,:),'Sample Name')),3,'first');
    Tcoli = find(~cellfun(@isempty,strfind(rawdata(1,:),'Target Name')),3,'first');
    Ccoli = find(~cellfun(@isempty,strfind(rawdata(1,:),'C')),3,'first');

    samples = unique(rawdata(2:end,Scoli));
    targets = unique(rawdata(2:end,Tcoli));

    Nsamples = numel(samples);
    Ntargets = numel(targets);

    % collect CT values

    CTmean = zeros([Nsamples Ntargets]);
    CTstd = zeros([Nsamples Ntargets]);

    for si = 1:Nsamples
        for ti = 1:Ntargets

            ind = strcmp(rawdata(:,Scoli), samples{si}) & strcmp(rawdata(:,Tcoli), targets{ti});
            if any(ind)
                %disp([samples{si} '; ' targets{ti}]);
                CTmean(si,ti) = str2double(rawdata{find(ind,1),Ccoli(2)});
                CTstd(si,ti) = str2double(rawdata{find(ind,1),Ccoli(3)});
            end
        end
    end

    data = struct('rawdata',{rawdata},...
        'targets',{targets},'samples',{samples},'CT',CTmean,'CTstd',CTstd,...
        'Nsamples', Nsamples, 'Ntargets', Ntargets);
end
