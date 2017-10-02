function upperleft = stitchedPreviews(dataDir, meta, type)

if ~exist('type', 'var')
	type = 'MIP';
end
    
MIPfiles = dir(fullfile(dataDir,'MIP','*MIP_*tif'));
s = strsplit(MIPfiles(1).name,['_' type]);
barefname = s{1};

gridSize = meta.montageGridSize;
pixelOverlap = round(meta.xSize*meta.montageOverlap/100);
posPerCondition = meta.posPerCondition;

% subsampling factor
if meta.xSize <= 1024
    ss = 2;
else
    ss = 4;
end

upperleft = {};

for wellnr = 13:meta.nWells
    
    conditionPositions = posPerCondition*(wellnr-1)+1:posPerCondition*wellnr;
    if isempty(gridSize) 
        gridSize = [round(posPerCondition/2) 2];
    end
    
    %nChannels = numel(MIPfiles)/(meta.nWells*posPerCondition);
    nChannels = meta.nChannels;
    
    for ci = 1:nChannels
        
        disp(['-------------processing channel ' num2str(ci) '--------']);
        
        imgs = {};
        tmax = meta.nTime;
        
        for ti = 1:tmax
            
            disp(['processing time ' num2str(ti)]);

            for pi = conditionPositions

                disp(['reading MIP ' num2str(pi)]);
                % gridSize 1 and 2 may be swapped, I have no way of knowing right now
                [i,j] = ind2sub(gridSize, pi - conditionPositions(1) + 1);

                fname = fullfile(dataDir,'MIP',[barefname sprintf(['_' type '_p%.4d_w%.4d.tif'],pi-1,ci-1)]);
                imgs{j,i} = imread(fname,ti);
            end

            % stitch together
            if ci == 1
                if ~isempty(pixelOverlap) %% && ti == 1
                    % get register positions of upper left corner
                    upperleft{wellnr}{ti} = registerImageGrid(imgs, pixelOverlap);
                else %if ti == 1 && isempty(pixelOverlap)
                    for pi = conditionPositions
                        [i,j] = ind2sub(gridSize,pi - conditionPositions(1) + 1);
                        upperleft{wellnr}{ti}{j,i} = [1+(j-1)*(meta.ySize + 50), 1+(i-1)*(meta.xSize + 50)];
                    end
                end
            end
            [stitched, upperleft{wellnr}{ti}] = stitchImageGrid(upperleft{wellnr}{ti}, imgs);
            % stitchImageGrid shift upperleft so all images are completely
            % within the stitched image 
            % CAUTION : upperleft = [y x]

            % make clean preview (not for quantitative analysis)
            small = imfilter(stitched,ones(ss)/ss^2);
            small = small(1:ss:end,1:ss:end);
            if ti == 1
                Ilim = stretchlim(small);
                Imin = double(min(small(small>0)));
                Imax = round(Ilim(2)*(2^16-1));
            end
            small = mat2gray(small, [Imin Imax]);
            small = uint16((2^16-1)*small);

            if ti == 1
                preview = zeros([size(small) tmax],'uint16');
            end
            preview(1:size(small,1), 1:size(small,2), ti) = small;
        end

        fname = fullfile(dataDir, [sprintf('stichedPreview_w%.4d_well',ci) num2str(wellnr) '.tif']);
        imwrite(preview(:,:,1), fname);
        for ti = 2:tmax
            imwrite(preview(:,:,ti), fname,'WriteMode','Append');
        end
    end
end
end