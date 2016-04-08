function displayPositions_Andor(meta)

    XYZ = meta.XYZ;

    %scatter(XYZ(:,1), XYZ(:,2))

    for i = 1:size(XYZ,1)
        text(XYZ(i,1), XYZ(i,2),num2str(i))
        w = 1024*meta.xres;
        h = 1024*meta.yres;
        rectangle('Position',[XYZ(i,1)-w/2,XYZ(i,2)-h/2,w,h])
    end
    axis([min(XYZ(:,1))-w max(XYZ(:,1))+w min(XYZ(:,2))-h max(XYZ(:,2))+h])
    axis equal
    axis off

    % XYZmean = mean(XYZ);
    % hold on
    % scatter(XYZmean(:,1), XYZmean(:,2),'r')
    % hold off
end