function stitchedPreviewsM(dataDirs, meta, type)

% make stitched previews from multiple recordings 
% (e.g. when the microscope was stopped to refocus)

if ~exist('type', 'var')
	type = 'MIP';
end

gridSize = meta.montageGridSize;
pixelOverlap = round(meta.xSize*meta.montageOverlap/100);
posPerCondition = meta.posPerCondition;

% subsampling factor
if meta.xSize <= 1024
    ss = 2;
else
    ss = 4;
end

for wellnr = 1:meta.nWells
    
    conditionPositions = posPerCondition*(wellnr-1)+1:posPerCondition*wellnr;
    if isempty(gridSize) 
        gridSize = [round(posPerCondition/2) 2];
    end
    
    nChannels = meta.nChannels;
    %nChannels = numel(MIPfiles)/(meta.nWells*posPerCondition);
    
    for ci = 1:nChannels
        
        disp(['-------------processing channel ' num2str(ci) '--------']);
        
        imgs = {};
        
        for di = 1:numel(dataDirs)
            
            MIPfiles = dir(fullfile(dataDirs{di},'MIP','*MIP_*tif'));
            s = strsplit(MIPfiles(1).name,['_' type]);
            barefname = s{1};

            fname = fullfile(dataDirs{di},'MIP',[barefname sprintf(['_' type '_p%.4d_w%.4d.tif'],0,0)]);
            tmax = numel(imfinfo(fname));
        
            for ti = 1:tmax

                disp(['processing time ' num2str(ti)]);

                for pi = conditionPositions

                    disp(['reading MIP ' num2str(pi)]);
                    % gridSize 1 and 2 may be swapped, I have no way of knowing right now
                    [i,j] = ind2sub(gridSize, pi - conditionPositions(1) + 1);

                    fname = fullfile(dataDirs{di},'MIP',[barefname sprintf(['_' type '_p%.4d_w%.4d.tif'],pi-1,ci-1)]);
                    imgs{j,i} = double(imread(fname,ti));
                end

                % stitch together
                if ci == 1
                    if ti == 1 && ~isempty(pixelOverlap)
                        % get register positions of upper left corner
                        upperleft = registerImageGrid(imgs, pixelOverlap);
                    elseif ti == 1 && isempty(pixelOverlap)
                        upperleft = {};
                        for pi = conditionPositions
                            [i,j] = ind2sub(gridSize,pi - conditionPositions(1) + 1);
                            upperleft{j,i} = [1+(j-1)*(meta.ySize + 50), 1+(i-1)*(meta.xSize + 50)];
                        end
                    end
                end
                stitched = stitchImageGrid(upperleft, imgs);

                % make clean preview (not for quantitative analysis)
                small = imfilter(stitched,ones(ss)/ss^2);
                small = small(1:ss:end,1:ss:end);
                if ti == 1 && di == 1
                    Ilim = stretchlim(small);
                    Imin = double(min(small(small>0)));
                    Imax = round(Ilim(2)*(2^16-1));
                end
                small = mat2gray(small, [Imin Imax]);
                small = uint16((2^16-1)*small);

                if ti == 1
                    preview = zeros([size(small) tmax],'uint16');
                end
                preview(:,:,ti) = small;
            end
            
            previews{di} = preview;
        end
        
        fname = fullfile(dataDirs{1}, [sprintf('stichedPreview_w%.4d_well',ci) num2str(wellnr) '.tif']);
        
        for di = 1:numel(dataDirs)
            
            MIPfiles = dir(fullfile(dataDirs{di},'MIP','*MIP_*tif'));
            s = strsplit(MIPfiles(1).name,['_' type]);
            barefname = s{1};

            fnametmp = fullfile(dataDirs{di},'MIP',[barefname sprintf(['_' type '_p%.4d_w%.4d.tif'],0,0)]);
            tmax = numel(imfinfo(fnametmp));
            
            preview = previews{di};
            if di == 1
                imwrite(preview(:,:,1), fname);
                tmin = 2;
            else
                tmin = 1;
            end
            for ti = tmin:tmax
                imwrite(preview(:,:,ti), fname,'WriteMode','Append');
            end
        end
    end
end
end