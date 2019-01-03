function colonies = dirtySegLive(dataParam, findColParam, segMaskParam, colonies)

for coli = dataParam.colonies 

    % load data
    fname = sprintf(dataParam.fnameformat, coli);
    img = readStack(fullfile(dataParam.dataDir,fname));
    
    % load segmentation
    segfname = sprintf(dataParam.segfnameformat, coli);
    seg = permute(squeeze(h5read(fullfile(dataParam.dataDir,segfname), '/exported_data')),[2 1 3]);
    seg = seg == 2;
    
    % make radial profiles over time
    firstfound = false;
    for ti = 1:dataParam.tmax

        fprintf('.');
        if mod(ti,80)==0
            fprintf('\n');
        end

        mask = seg(:,:,ti);
        
        range = [];
        [coloniesTmp, cleanmask, welllabel] = findColonies(mask, range, dataParam.meta, findColParam);

        if ti == 1
            colonies(coli) = coloniesTmp;
        else
            colonies(coli).boundingBox = cat(1, colonies(coli).boundingBox, coloniesTmp.boundingBox);
            colonies(coli).center = cat(1, colonies(coli).center, coloniesTmp.center);
        end

        % sometimes the first frame is a fail so you can't populate the
        % structure with that
        if ~firstfound && ~isempty(coloniesTmp.radius) 
            colonies(coli).radiusPixel = coloniesTmp.radiusPixel;
            colonies(coli).radiusMicron = coloniesTmp.radiusMicron;
            firstfound = true;
        end
        
    	% save colony mask
        if ti == 1
            imwrite(cleanmask, sprintf(fullfile(dataParam.dataDir,'cleanmask_col%d.tif'),coli),'Compression','none');
        else
            imwrite(cleanmask, sprintf(fullfile(dataParam.dataDir,'cleanmask_col%d.tif'),coli),'Compression','none','WriteMode','append');
        end

        % radial profile
        if ~isnan(colonies(coli).boundingBox(end,1))

            % make masks
            if strcmp(dataParam.dataType, 'SMAD4')                
                [N, D] = dirtySegMasksSMAD4(cleanmask, mask, segMaskParam);
            elseif strcmp(dataType, 'BCAT') 
                [N, D] = dirtySegMasksBCAT(cleanmask, mask, segMaskParam);
            else
                error('unknown data type');
            end
            
            % dealing with clipped colonies
            b = colonies(coli).boundingBox(ti,:);
            b(3) = max(b(3),1); b(1) = max(b(1),1);
            b(4) = min(b(4), size(img,1)); b(2) = min(b(2), size(img,2));
            N = N(b(3):b(4),b(1):b(2));
            D = D(b(3):b(4),b(1):b(2));
            colimg = img(b(3):b(4),b(1):b(2),:,:,ti);

            % save diagnostic image
            if ti == 1 || ti == tmax
                I = imadjust(mat2gray(colimg));
                s = 0.2;
                RGB = cat(3, I + s*N, I + s*D ,I);
                imwrite(RGB, sprintf(fullfile(dataParam.dataDir,'segmentation_col%d_t%d.tif'),coli,ti),'Compression','none','WriteMode','append');
            end
            
            % make radial profile
            colonies(coli).makeRadialAvgNoSeg(colimg, N, D, dataParam.meta.colMargin, ti)
        else
            % create and empty struct to keep track of missed time pts
            if ti>1
                colonies(coli).radialProfile(ti) = struct('BinEdges',[],'NucAvg',[],'NucStd',[],'CytAvg',[],'CytStd',[]);
            else
                colonies(coli).radialProfile = struct('BinEdges',[],'NucAvg',[],'NucStd',[],'CytAvg',[],'CytStd',[]);
            end
        end
    end
toc
end

end