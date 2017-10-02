function [cytstart, cytend] = getCyt(results, data)

    frapframe = results.frapframe;
    imstart = squeeze(data(:,:,1));
    imend = squeezed(data(:,:,min(results.tmax)));
    imfrap = squeeze(data(:,:,frapframe));

    cytmask = {};
    bg = {};
    nuc = {};
    cyt = {};
    
    im = imstart;
    imshow(imadjust(im),[])
    hold on
    for shapeIdx = 1:numel(x)
        cytstart{shapeIdx} = mean(im(cytmask{shapeIdx}));
    end
    hold off
    
    im = imend;
    imshow(imadjust(im),[])
    hold on
    for shapeIdx = 1:numel(x)
        
        plot(x{shapeIdx}, y{shapeIdx});
        xp = x{shapeIdx} - 0.3*(x{shapeIdx}-nanmean(x{shapeIdx}));
        yp = y{shapeIdx} - 0.3*(y{shapeIdx}-nanmean(y{shapeIdx}));
        
        nucmask{shapeIdx} = poly2mask(xp, yp, size(im,2), size(im,1));
        %cytmask{shapeIdx} = roipoly;
        
        bg{shapeIdx} = mean(imfrap(nucmask{shapeIdx}));
        nucend{shapeIdx} = mean(im(nucmask{shapeIdx}));
        cytend{shapeIdx} = mean(im(cytmask{shapeIdx}));
    end
    hold off
    
end