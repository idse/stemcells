function [traces, tracesnorm, maskall] = readFRAPprofile(data, x,y)

    tcut = size(data,3);
    
    traces = zeros([numel(x) tcut]);
    tracesnorm = zeros([numel(x) tcut]);
    maskall = false([size(data,1) size(data,2)]);

    for shapeIdx = 1:numel(x)

        mask = poly2mask(x{shapeIdx}, y{shapeIdx}, size(data,2), size(data,1));
        
        pixelArea = sum(mask(:));
        s = round(sqrt(pixelArea)/8);
        
        mask = imerode(mask, strel('disk', s));
        maskall = maskall | mask;

        val = zeros([1 tcut]);
        for ti = 1:tcut
            im = data(:,:,ti);
            val(ti) = mean(im(mask));
        end

        traces(shapeIdx,:) = val;
        tracesnorm(shapeIdx,:) = (val - min(val))/(max(val) - min(val));
    end
end