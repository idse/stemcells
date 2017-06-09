function [x,y, shapeTypes] = readFRAPregions(omeMeta)

    ROIidx = 0; % always only one ROI
    x = {};
    y = {};
    shapeTypes = {};
    Nfrapped = omeMeta.getShapeCount(ROIidx)/2;

    for shapeIdx = 1:Nfrapped

        shapeType = char(omeMeta.getShapeType(ROIidx,shapeIdx-1));
        shapeTypes{shapeIdx} = shapeType;

        if strcmp(shapeType, 'Polygon')

            disp(['reading polygon ' num2str(shapeIdx)]);

            s = strsplit(char(omeMeta.getPolygonPoints(ROIidx,shapeIdx-1)),{' ',','});
            s = cellfun(@str2double, s,'UniformOutput',false);
            s = cat(1, s{:});

            x{shapeIdx} = s(1:2:end);
            y{shapeIdx} = s(2:2:end);

        elseif strcmp(shapeType,'Ellipse')

            disp('ellipse?');
            X = omeMeta.getEllipseX(ROIidx, shapeIdx-1);
            Y = omeMeta.getEllipseY(ROIidx, shapeIdx-1);
            RX = omeMeta.getEllipseRadiusX(ROIidx, shapeIdx-1);
            RY = omeMeta.getEllipseRadiusY(ROIidx, shapeIdx-1);
            A = omeMeta.getEllipseTransform(ROIidx, shapeIdx-1);
            % this is a pain

        elseif strcmp(shapeType, 'Rectangle')

            rx = double(omeMeta.getRectangleX(ROIidx, shapeIdx-1));
            ry = double(omeMeta.getRectangleY(ROIidx, shapeIdx-1));
            rh = double(omeMeta.getRectangleHeight(ROIidx, shapeIdx-1));
            rw = double(omeMeta.getRectangleWidth(ROIidx, shapeIdx-1));

            x{shapeIdx} = [rx, rx + rw, rx + rw, rx, rx];
            y{shapeIdx} = [ry, ry, ry + rh, ry + rh, ry];
        end
    end
end