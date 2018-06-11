function visualizePressure(g, options)
    
    L = options.L;

    options.edgeColor = sumTension(g);  %,'colorTable',g.A);
    options.transparent = false;
    
    CT = g.stress;
    CT(isnan(CT)) = 0;
    options.colorTable = CT;

    cells = g.cells;
    bonds = g.bonds;
    verts = g.verts;

    nCells = numel(cells);
    maxNN = 14;
    faces = nan([nCells maxNN]);
    for ci = 1:length(cells) 
        C = cells{ci};
        faces(ci,1:numel(C)) = bonds(C,1);
    end
    
    minVal = options.pDispRange(1);
    maxVal = options.pDispRange(2);
    
    cmap = jet(256);
    colorsIdx = round(255*min(CT - minVal, maxVal-minVal)/(maxVal-minVal) + 1);
    colorsIdx(colorsIdx < 1) = 1;
    colorsIdx(colorsIdx > 256) = 256;
    faceCol = cmap(colorsIdx,:);
    
    clf
    patch('Faces',faces,'Vertices',verts,'FaceVertexCData',faceCol);
    shading faceted
    axis equal
    axis off

    % clone outline
    hold on
    %cloneCells = find(~isnan(g.clones));
    cloneBonds = [g.cells{~isnan(g.clones)}]';
    cloneEdge = cloneBonds(isnan(g.clones(g.bonds(cloneBonds,4))));
    V1 = g.verts(g.bonds(cloneEdge,1),:);
    V2 = g.verts(g.bonds(cloneEdge,2),:);
    for i = 1:size(V1,1)
        plot([V1(i,1) V2(i,1)],[V1(i,2) V2(i,2)],'k','LineWidth',3);
    end
	hold off          

    axis equal
    
    axis([-L L -L L])
    axis off
    set(gcf,'Color','white');
end