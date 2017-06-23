function [colonies, cleanmask, welllabel] = findColonies(mask, range, meta, s)
    % find colonies in binary mask
    %
    % [colonies, cleanmask] = findColonies(mask, range, meta, s)
	%
    % mask:         binary image (thresholded DAPI)
    % range:        [xmin xmax ymin ymax] into mask
    % meta:         meta data object
    % s:            radius of strel for cleaning up, in micron
    %
    % colonies:     array of colony objects
    % cleanmask:    cleaned up binary mask
    
    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------
    
    colRadiusPixel = meta.colRadiiPixel;
    
    % clean up mask
    %-----------------------

    minArea = floor(pi*min(colRadiusPixel).^2/2);
    maxArea = ceil(2*pi*max(colRadiusPixel).^2);

    % try to make solid objects out of colonies
    cleanmask = mask(range(3):range(4),range(1):range(2));
    
    cleanmask = imclose(cleanmask,strel('disk',s));
    cleanmask = imfill(cleanmask,'holes');
    cleanmask = imopen(cleanmask,strel('disk',s));
    cleanmask = bwareaopen(cleanmask,minArea);
    
    %remove the colonies to get the mask of the wells
    wellmask = ~bwareaopen(cleanmask,maxArea);
    
    % remove colonies that are not round
    CC = bwconncomp(cleanmask);
    stats = regionprops(CC, 'Eccentricity', 'BoundingBox','Centroid');
    cleanmask(cat(1,CC.PixelIdxList{[stats.Eccentricity] > 0.5})) = false;
    
    %labelled well image
    welllabel = bwlabel(wellmask);
    
    
    % find colonies
    %-----------------------
    
    goodColIdx = [stats.Eccentricity] < 0.5;
    nColonies = sum(goodColIdx);
    if nColonies == 0
        error('no complete colonies found');
    end
    
    % determine colony size
    
 
    bb = cat(1,stats.BoundingBox);
    bb = bb(goodColIdx,:);
    linsize = max(bb(:,3:4),[],2);
    
    if numel(colRadiusPixel) > 1
    edges = imfilter(2*colRadiusPixel,[1 1]/2);
    edges = [0, edges(1:end-1), 2*edges(end-1)-edges(end-2)];
    [~,colType] = histc(linsize,edges);
    else
        colType = ones(size(linsize));
    end
    
    % colony centers
    CM = cat(1,stats.Centroid);
    CM = CM([stats.Eccentricity] < 0.5,:);
    
    %get wells
     well = welllabel(sub2ind(size(welllabel),floor(CM(:,2)),floor(CM(:,1))));

    
    % shift to absolute position
    CM(:,1) = CM(:,1) + double(range(1) - 1);
    CM(:,2) = CM(:,2) + double(range(3) - 1);

    
    % note: colRadii here is the radii of individual colonies
    % meta.colRadii contains the small number of possible radii
    
    colRadii = zeros(size(colType)); colRadiiMicron=colRadii;
    colRadii(colType > 0) = colRadiusPixel(colType(colType>0));
    colRadiiMicron(colType>0) = meta.colRadiiMicron(colType(colType>0))
    % sort by radius
    [colRadii, idx] = sort(cat(1,colRadii),'descend');
    CM = CM(idx,:);
    colRadiiMicron = colRadiiMicron(idx);
    
    % colony ranges (improved bounding box of identical size)
    colMargin = 10;
    Rmax = ceil(colRadii) + colMargin;
    colxmin = floor(CM(:,1)) - Rmax;
    colxmax = ceil(colxmin) + 2*Rmax;
    colymin = floor(CM(:,2)) - Rmax;
    colymax = ceil(colymin) + 2*Rmax;
    colrange = [colxmin colxmax colymin colymax];

    % make that the colonies are completely within range
    contained = (colxmin > range(1)) & (colxmax < range(2)) &...
                               (colymin > range(3)) & (colymax < range(4));
    CM = CM(contained,:);
    colRadii = colRadii(contained);
    colrange = colrange(contained,:);
    colRadiiMicron = colRadiiMicron(contained);
    nColonies = sum(contained);
                           
    % array of colony objects
    colonies(nColonies) = Colony;
    for i = 1:nColonies;
        % Colony(nChannels, center, radiusPixel, radiusMicron, boundingBox)
        colonies(i) = Colony(meta.nChannels, CM(i,:), colRadii(i),...
                                        colRadiiMicron(i), colrange(i,:),well(i)); 
    end
end