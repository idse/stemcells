function makeStitchedStack(dataDir, meta, upperleft, wells, margin)

    if ~exist('wells','var')
        wells = 1:meta.nWells;
    end
    if ~exist('margin','var')
        margin = 0;
    end
    if strfind(meta.filename, '_z')
        zsplit = true;
    else
        zsplit = false;
    end
    if strfind(meta.filename, '_w')
        channelsplit = true;
    else
        channelsplit = false;
    end
    
    gridSize = meta.montageGridSize;
    
    % order to save should be xyczt
    % so it should be possible to not fill up the memory by only reading
    % the required timepoints
    
    if zsplit && ~channelsplit
        
        for ti = 1:meta.nTime
        
            imgs = {};
            
            for wellnr = wells

                if ti == 1
                    disp(['processing well/colony: ' num2str(wellnr)]);
                else
                    fprintf('.');
                end
            
                a = (wellnr-1)*meta.posPerCondition;
                conditionPositions = a:a+meta.posPerCondition-1;

                for zi = 0:meta.nZslices-1
                    
                    for ci = 1:meta.nChannels
                        
                        outfname = sprintf('col_id%d_c%d.tif', wellnr, ci);
                        
                        for pi = conditionPositions

                            [i,j] = ind2sub(gridSize,pi - conditionPositions(1) + 1);
                            fname = fullfile(dataDir, sprintf(meta.filename, pi, zi));

                            % xyct
                            frameidx = meta.nChannels*(ti - 1) + ci;
                            imgs{j,i,ci} = imread(fname, frameidx);
                        end
                        
                        %disp(['ti ' num2str(ti) ', zi ' num2str(zi) ', ci ' num2str(ci)]);
                        stitched = stitchImageGrid(upperleft{wellnr}{ti}, imgs(:,:, ci));

                        % this only works if each frame is the same size
                        if ti==1 && zi==0 && ci == 1
                            stitchedSize = size(stitched);
                        end
                        
                        ymin = margin + 1;
                        ymax = min(size(stitched,1),stitchedSize(1)-margin);
                        
                        xmin = margin + 1;
                        xmax = min(size(stitched,2),stitchedSize(2)-margin);
                        
                        stitchedp = zeros(stitchedSize-2*margin, 'uint16');
                        stitchedp(1:ymax-margin,1:xmax-margin) = stitched(ymin:ymax,xmin:xmax);

                        if ti==1 && zi==0
                            imwrite(stitchedp, fullfile(dataDir,outfname));
                        else
                            imwrite(stitchedp, fullfile(dataDir,outfname), 'WriteMode','Append');
                        end
                    end
                end
            end
        end
        
    else
        
        error('implement this case');
    end
end