function results = getCyt(results, data)

    cytstart = zeros([numel(results.x) size(data,3)]);    
    cytend = zeros([numel(results.x) size(data,3)]); 
    
    cytmask = {};
    imshow(imadjust(data(:,:,1)),[])
    hold on
    for shapeIdx = 1:numel(results.x)

        plot(results.x{shapeIdx},results.y{shapeIdx});
        cytmask{shapeIdx} = roipoly;
        
        
        for ti = 1:size(data,3)
            im = data(:,:,ti);
            cytstart(shapeIdx,ti) = mean(im(cytmask{shapeIdx}));
        end
    end
    hold off
    
    cytmask = {};
    imshow(imadjust(data(:,:,min(results.tmax))),[])
    hold on
    for shapeIdx = 1:numel(results.x)
        
        plot(results.x{shapeIdx},results.y{shapeIdx});
        cytmask{shapeIdx} = roipoly;
        for ti = 1:size(data,3)
            im = data(:,:,ti);
            cytend(shapeIdx,ti) = mean(im(cytmask{shapeIdx}));
        end
    end
    hold off 
    close;
    
    results.cytstart = cytstart;
    results.cytend = cytend;
end