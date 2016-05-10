function stitched = stitchImageGrid(upperleft, imgs)
    % combine image grid into single image
    %
    % stitched = stitchImageGrid(upperleft, imgs)
    %
    %
    % stitched:     combined image
    %
    % upperleft:    cell array of positions of upperleft corners with
    %               upperleft corner of upperleft image being (1,1)
    %               as produced by registerImageGrid(..)
    % imgs:         cell array of images 
    %               with rows and cols corresponding to grid
    %
    % see also registerImageGrid

    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------

    % note: could be made fancier by averaging the overlap

    UL = cat(1,upperleft{:});
    dUL = 1 - min(UL);
    for i = 1:numel(upperleft)
        if ~isempty(upperleft{i})
            upperleft{i} = upperleft{i} + dUL;
        end
    end

    % again assuming square images
    N = size(imgs{1,1},1);
    totalSize = max(UL) + N;

    stitched = zeros(totalSize,'uint16');
    for i = 1:size(imgs,1)
        for j = 1:size(imgs,2)
            if ~isempty(upperleft{i,j})
                I = upperleft{i,j}(1):upperleft{i,j}(1)+N-1;
                J = upperleft{i,j}(2):upperleft{i,j}(2)+N-1;
                stitched(I,J) = imgs{i,j};
            end
        end
    end
end