function res = readFRAPprofile2(data, omeMeta, res)

    % res = results structure
    tcut = size(data,3);
    
    if ~isfield(res, 'cytx')
        res.nucxstart = {};
        res.nucxend = {};
        res.cytxstart = {};
        res.cytxend = {};
        res.nucx = {};
        res.cytx = {};
        res.shrink = 0.5; 
    end

    % read FRAP regions
    [x,y,shapeTypes] = readFRAPregions(omeMeta);
    results.shapeTypes = shapeTypes;

    % get all the masks
    for shapeIdx = 1:numel(x)
        
        if isempty(res.nucxend) || numel(res.nucxend) < shapeIdx
        
            disp('select final nuclear mask');
            imshow(imadjust(data(:,:,end-1)),[],'InitialMagnification',200)
            hold on
            plot(x{shapeIdx},y{shapeIdx},'LineWidth',2);
            h = impoly(gca);
            res.nucxstart{shapeIdx} = [x{shapeIdx},y{shapeIdx}];
            res.nucxend{shapeIdx} = getPosition(h);
            hold off

            disp('select initial cytoplasmic mask');
            imshow(imadjust(data(:,:,1)),[])
            hold on
            plot(x{shapeIdx},y{shapeIdx},'LineWidth',2);
            h = impoly(gca);
            res.cytxstart{shapeIdx} = getPosition(h);
            hold off

            disp('select final cytoplasmic mask');
            imshow(imadjust(data(:,:,end-1)),[])
            hold on
            plot(x{shapeIdx},y{shapeIdx},'LineWidth',2);
            plot(res.cytxstart{shapeIdx}(:,1), res.cytxstart{shapeIdx}(:,2));
            h = impoly(gca);
            res.cytxend{shapeIdx} = getPosition(h);
            hold off
            
            close;
        end
        
        % dynamic polygons
        res.nucx{shapeIdx} = polymorph(res.nucxstart{shapeIdx}, res.nucxend{shapeIdx}, tcut);
        res.cytx{shapeIdx} = polymorph(res.cytxstart{shapeIdx}, res.cytxend{shapeIdx}, tcut);
    end

    % read the profile
    tracesNuc = zeros([numel(x) tcut]);
    tracesNucNorm = zeros([numel(x) tcut]);
    tracesCyt = zeros([numel(x) tcut]);
    tracesCytNorm = zeros([numel(x) tcut]);
    
    for shapeIdx = 1:numel(x)
        
        nucval = zeros([1 tcut]);
        cytval = zeros([1 tcut]);
        
        for ti = 1:tcut

            nucp = res.nucx{shapeIdx}(:,:,ti);
            nucp = nucp - res.shrink*(nucp - mean(nucp));
            
            cytp = res.cytx{shapeIdx}(:,:,ti);
            %cytp = cytp - prof.shrink*(cytp - mean(cytp));
            
            nucmask = poly2mask(nucp(:,1), nucp(:,2), size(data,2), size(data,1));
            cytmask = poly2mask(cytp(:,1), cytp(:,2), size(data,2), size(data,1));
            
            im = data(:,:,ti);
            nucval(ti) = mean(im(nucmask));
            cytval(ti) = mean(im(cytmask));
        end
        
        bg = min(nucval);
        
        tracesNuc(shapeIdx,:) = nucval;
        tracesNucNorm(shapeIdx,:) = (nucval - bg)/(mean(nucval(1:2)) - bg);
        
        tracesCyt(shapeIdx,:) = cytval;
        tracesCytNorm(shapeIdx,:) = (cytval - bg)/(mean(cytval(1:2)) - bg);
    end

    res.tracesNuc = tracesNuc;
    res.tracesNucNorm = tracesNucNorm;
    res.tracesCyt = tracesCyt;
    res.tracesCytNorm = tracesCytNorm;
end