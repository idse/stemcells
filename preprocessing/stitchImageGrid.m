function stitched = stitchImageGrid(upperleft, imgs)

    % note: could be made fancier by averaging the overlap

    UL = cat(1,upperleft{:});
    dUL = 1 - min(UL);
    for i = 1:numel(upperleft)
        upperleft{i} = upperleft{i} + dUL;
    end

    % again assuming square images
    N = size(imgs{1,1},1);
    totalSize = max(UL) + N;

    stitched = zeros(totalSize,'uint16');
    for i = 1:4
        for j = 1:4
            I = upperleft{i,j}(1):upperleft{i,j}(1)+N-1;
            J = upperleft{i,j}(2):upperleft{i,j}(2)+N-1;
            stitched(I,J) = imgs{i,j};
        end
    end
end