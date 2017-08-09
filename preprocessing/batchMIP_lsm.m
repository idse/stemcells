function batchMIP_lsm(inputdir, outputdir, channels, saveidx, tmax, positions)

    dirstruct = cat(1,  dir(fullfile(inputdir,'*oib')),...
                        dir(fullfile(inputdir,'*oif'))  );
                    
    if ~isempty(dirstruct)
        % just a bunch of oib file case (FRAP)
        
        if ~exist(outputdir,'dir')
            mkdir(outputdir);
        end

        oibfiles = {dirstruct.name};
    
        for fi = 1:numel(oibfiles)
            
            oibfile = oibfiles{fi};
            [~,barefname,~] = fileparts(oibfile);
            barefname = strrep(barefname,'.','dot');

            % read the data
            [img, ~] = readStack2(fullfile(inputdir,oibfile));

            % save nuclear MIP
            for ci = 1:numel(channels)
                MIP = uint16(squeeze(sum(img(:,:,channels(ci)+1,:,:),4)));
                fname = fullfile(outputdir, sprintf([barefname '_MIP_w%.4d.tif'], channels(ci)));
                if exist(fname,'file')
                    delete(fname);
                end
                imwrite(MIP(:,:,1), fname);
                for i = 2:size(MIP,3)
                    imwrite(MIP(:,:,i), fname,'WriteMode','Append');
                end
            end
        end

    else
        % multiposition oif case
        
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
end