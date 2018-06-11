function visualizeImageStitchQuality(imgsNuc, pixelOverlap)
    % future versions of registerImageGrid should stitch together images
    % taking into account which overlaps actually contain enough
    % information to determine a faitful displacement vector, this function
    % visualizes this as a starting point for improvements

    gridUpperleft = registerImageGrid(imgsNuc, []);
    [~,links] = registerImageGrid(imgsNuc, pixelOverlap);
    gridImg = stitchImageGrid(gridUpperleft, imgsNuc);
    gridImg = imadjust(mat2gray(gridImg));
    imshow(gridImg)
    N = 1024 + 50;
    for i = 1:size(imgsNuc,1)
        for j = 1:size(imgsNuc,2)-1

            % horizontal links
            x = 3*N/4 + (j-1)*N;
            y = N/2 + (i-1)*N;

            color = 'green';
            if links(2*i-1,2*j) < 1.05
                color = 'red';
            end
            line([x x+N/2],[y y],'LineWidth',2,'Color',color)
        end
    end
    for i = 1:size(imgsNuc,1)-1
        for j = 1:size(imgsNuc,2)
            % vertical links
            x = N/2 + (j-1)*N;
            y = 3*N/4 + (i-1)*N;

            color = 'green';
            if links(2*i,2*j-1) < 1.05
                color = 'red';
            end
            line([x x],[y y+N/2],'LineWidth',2,'Color',color)
        end
    end
end