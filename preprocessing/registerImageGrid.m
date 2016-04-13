function upperleft = registerImageGrid(imgs, pixelOverlap)
    % register a grid of overlapping images
    % 
    % upperleft = registerImageGrid(imgs, pixelOverlap)
    %
    %
    % upperleft:    cell array of positions of upperleft corners with
    %               upperleft corner of upperleft image being (1,1)
    %
    % imgs:         cell array of images 
    %               with rows and cols corresponding to grid
    % pixelOverlap: width of overlapping strip in pixels

    % note: could be made fancier by combining redundant shift information
    
    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------
    
    % I will assume NxN square images
    N = size(imgs{1,1},1);

    belowshift = {};
    for i = 1:size(imgs, 1)-1     
        for j = 1:size(imgs, 2)
            img = imgs{i,j}(end-pixelOverlap+1:end,:);
            below = imgs{i+1,j}(1:pixelOverlap,:);
            [shifti,shiftj] = xcorr2fft(below, img);
            belowshift{i+1,j} = [shifti,shiftj];
        end
    end

    rightshift = {};
    for j = 1:size(imgs, 2)-1
        for i = 1:size(imgs, 1)
            img = imgs{i,j}(:,end-pixelOverlap+1:end);
            right = imgs{i,j+1}(:,1:pixelOverlap);
            [shifti,shiftj] = xcorr2fft(right, img);
            rightshift{i,j+1} = [shifti,shiftj];
        end
    end

    Np = N - pixelOverlap;

    upperleft = {};
    upperleft{1,1} = [1 1];

    % register each row in the first column
    j = 1;
    for i = 2:size(imgs, 1)
        shift = belowshift{i,j};
        upperleft{i,j} = upperleft{i-1,j} + [Np + shift(1) + 1, shift(2)];
    end

    % for each row, register all the columns relative to the first one
    for i = 1:size(imgs, 1)
        for j = 2:size(imgs, 2)
            shift = rightshift{i,j};
            upperleft{i,j} = upperleft{i,j-1} + [shift(1), Np + shift(2) + 1];
        end
    end
end