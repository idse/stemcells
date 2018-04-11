function res = readFRAPprofileCytbleach(img, res, tmax)

    nucChannel = 2;
    S4Channel = 1;
    zi = 1;
    data = squeeze(img(:,:,S4Channel,zi,:));
    
    if ~exist('tmax','var') || isempty(tmax)
       tmax = size(data,3)-1; 
    end

    % res = results structure
    tcut = size(data,3);

    if ~isfield(res, 'cytx')
        res.nucxstart = {};
        res.nucxend = {};
        res.cytxstart = {};
        res.cytxend = {};
        res.nucx = {};
        res.cytx = {};
    end
    if ~isfield(res, 'shrink')
        res.shrink = 0.3; 
    end

    % get background levels
    Iraw = data(:,:,1);
    I = mat2gray(Iraw);

    mask = I > graythresh(I);
    mask = imclose(mask,strel('disk',3));
    mask = imopen(imfill(mask,'holes'), strel('disk',3));
    mask = bwareaopen(mask,1000);
    bgmask = imerode(~mask,strel('disk',40));
    if sum(bgmask(:)) < 100
        warning('cant determine background with standard margin');
        bgmask = imerode(~mask,strel('disk',10));
    end
    if sum(bgmask(:)) < 100
        warning('cant determine background');
        bgmask = imerode(~mask,strel('disk',10));
    end
    
    bg = mean(Iraw(bgmask));% std(double(Iraw(bgmask)))]
    res.bgempty = bg;
    
    R = imadjust(data(:,:,1));
    B = imadjust(data(:,:,res.frapframe));
    G = imadjust(squeeze(img(:,:,nucChannel,zi,res.frapframe)));
    initim = cat(3,R,R+G,R+B);
    
    R = imadjust(data(:,:,tmax));
    G = imadjust(squeeze(img(:,:,nucChannel,zi,tmax)));
    finim = cat(3,R,R+G,R);
        
    imshow(initim);
    Nfrapped = input('How many cells were bleached?\n');
    close;
    
    % get all the masks
    for shapeIdx = 1:Nfrapped

        disp('select initial nuclear mask');
        imshow(initim,'InitialMagnification',200)
        hold on
        h = impoly(gca);
        res.nucxstart{shapeIdx} = getPosition(h);
        x{shapeIdx} = res.nucxstart{shapeIdx}(:,1);
        y{shapeIdx} = res.nucxstart{shapeIdx}(:,2);
        hold off

        disp('select final nuclear mask');
        imshow(finim,'InitialMagnification',200)
        hold on
        plot(x{shapeIdx},y{shapeIdx},'LineWidth',2);
        h = impoly(gca);
        res.nucxend{shapeIdx} = getPosition(h);
        hold off

        disp('select initial cytoplasmic mask');
        imshow(initim)
        hold on
        plot(x{shapeIdx},y{shapeIdx},'LineWidth',2);
        h = impoly(gca);
        res.cytxstart{shapeIdx} = getPosition(h);
        hold off

        disp('select final cytoplasmic mask');
        imshow(finim)
        hold on
        plot(x{shapeIdx},y{shapeIdx},'LineWidth',2);
        plot(res.cytxstart{shapeIdx}(:,1), res.cytxstart{shapeIdx}(:,2));
        h = impoly(gca);
        res.cytxend{shapeIdx} = getPosition(h);
        hold off

        close;
        
        % dynamic polygons
        tic
        res.nucx{shapeIdx} = polymorph(res.nucxstart{shapeIdx}, res.nucxend{shapeIdx}, tcut);
        res.cytx{shapeIdx} = polymorph(res.cytxstart{shapeIdx}, res.cytxend{shapeIdx}, tcut);
        toc
    end

    % read the profile
    tracesNuc = zeros([Nfrapped tcut]);
    tracesNucNorm = zeros([Nfrapped tcut]);
    tracesCyt = zeros([Nfrapped tcut]);
    tracesCytNorm = zeros([Nfrapped tcut]);
    
    disp(['shrink mask :' num2str(res.shrink)]);
    
    for shapeIdx = 1:Nfrapped

        nucval = zeros([1 tcut]);
        cytval = zeros([1 tcut]);
        
        for ti = 1:tcut

            nucp = res.nucx{shapeIdx}(:,:,ti);
            nucp = nucp - res.shrink*(nucp - repmat(mean(nucp),[size(nucp,1) 1]));
            
            cytp = res.cytx{shapeIdx}(:,:,ti);
            %cytp = cytp - prof.shrink*(cytp - mean(cytp));
            
            nucmask = poly2mask(nucp(:,1), nucp(:,2), size(data,2), size(data,1));
            cytmask = poly2mask(cytp(:,1), cytp(:,2), size(data,2), size(data,1));
            
            im = data(:,:,ti);
            nucval(ti) = mean(im(nucmask));
            cytval(ti) = mean(im(cytmask));
        end
        
        %bg = min(nucval);
        
        tracesNuc(shapeIdx,:) = nucval;
        tracesCyt(shapeIdx,:) = cytval;
    end

    traceMin = min(min(tracesNuc));
    if traceMin < bg
        warning(['bg from empty space ' num2str(bg) ' was higher than FRAP curve ' num2str(traceMin)]);
        bg = traceMin;
    end

    % normalize the traces
    for shapeIdx = 1:Nfrapped
 
        nucval = tracesNuc(shapeIdx,:);
        tracesNucNorm(shapeIdx,:) = (nucval - bg)/(mean(nucval(1:2)) - bg);
        cytval = tracesCyt(shapeIdx,:);
        tracesCytNorm(shapeIdx,:) = (cytval - bg)/(mean(cytval(1:2)) - bg);
    end

    res.tracesNuc = tracesNuc;
    res.tracesNucNorm = tracesNucNorm;
    res.tracesCyt = tracesCyt;
    res.tracesCytNorm = tracesCytNorm;
    res.bg = bg;
end