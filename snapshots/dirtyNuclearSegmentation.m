function nucseg = dirtyNuclearSegmentation(im, opts)

    % parameters
    %-------------
    
    if isfield(opts, 's')
        s = opts.s;
    else
        s = 10;
        warning('filter scale was not provided, using default');
    end
    if isfield(opts, 'normalizedthresh')
        normalizedthresh = opts.normalizedthresh;
    else
        normalizedthresh = 0.5;
        warning('normalized threshold was not provided, using default');
    end
    if isfield(opts, 'absolutethresh')
        absolutethresh = opts.absolutethresh;
    else
        absolutethresh = graythresh(mat2gray(im))*max(im(:));%1000;
        warning('absolute threshold was not provided, using default');
    end
    if isfield(opts, 'areacutoff')
        areacutoff = opts.areacutoff;
    else
        areacutoff = s^2;
        warning('area cutoff was not provided, using default');
    end
    if ~isfield(opts, 'toobigNstd')
        opts.toobigNstd = 2;
    end

    % actual code
    %-------------
    
    if ~isfield(opts, 'mask')
        imdil = imdilate(im,strel('disk',2*s));
        imer = imerode(im,strel('disk',2*s));
        imnormalized = single(im - imer)./single(imdil - imer);
        clear imdil;
        clear imer;
    
        mask = imnormalized > normalizedthresh & im > absolutethresh;
        mask = imclose(mask, strel('disk',2));
        mask = bwareaopen(mask, areacutoff);
    else
        mask = opts.mask;
    end
    logim = imfilter(im, -fspecial('log',3*s, s));
    seeds = imextendedmax(logim, round(mean(logim(:))));
    
    % watershed
    D = bwdist(~mask);
    Dc = imcomplement(mat2gray(D));
    Ds = mat2gray(bwdist(seeds));
    Dtot = Ds + Dc;

    V = imimposemin(Dtot, seeds);
    L = watershed(V);
    L(~mask) = 0;

    %Lrgb = label2rgb(L,'jet','k','shuffle');
    %imshow(Lrgb)

    nucseg = imopen(L > 0,strel('disk',2));
    nucseg = bwareaopen(nucseg, areacutoff);
    
    %impp = imadjust(mat2gray(im));
    %newmaskedge = nucseg - imerode(nucseg,strel('disk',1));
    %overlay = cat(3, impp, newmaskedge, impp + seeds);
    %imshow(overlay)
    
    % remove largest nuclei (probably segmentation error)
    stats = regionprops(nucseg,'Area');
    Acutoff = round(mean([stats.Area]) + opts.toobigNstd*std([stats.Area]));
    disp(['A cutoff: ' num2str(Acutoff)]);
    toobig = bwareaopen(nucseg, Acutoff);
    nucseg(toobig) = 0;
end
