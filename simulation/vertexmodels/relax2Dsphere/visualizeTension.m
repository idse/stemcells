function visualizeTension(g, options)

    L = options.L;

    options.edgeColor = sumTension(g);  %,'colorTable',g.A);
    CT = g.clones;
    CT(isnan(CT)) = 0;
    options.colorTable = CT;
    options.transparentColor = 0;

    CT = options.edgeColor;
    minVal = options.edgeColorRange(1);
    maxVal = options.edgeColorRange(2);
    edgeCData = 64*min(CT - minVal, maxVal-minVal)/(maxVal-minVal);

    nCells = numel(g.cells);
    for ci = 1:nCells
 
        cdata = edgeCData(g.cells{ci});
        vertsi = g.bonds(g.cells{ci});
        patch(g.verts(vertsi,1), g.verts(vertsi,2), cdata,...
                        'EdgeColor','flat','FaceColor','None','LineWidth',2,...
                        'CDataMapping','direct');
    end

    % clone outline
    hold on
    cloneCells = find(~isnan(g.clones));
    cloneBonds = [g.cells{cloneCells}]';
    cloneEdge = cloneBonds(isnan(g.clones(g.bonds(cloneBonds,4))));
    V1 = g.verts(g.bonds(cloneEdge,1),:);
    V2 = g.verts(g.bonds(cloneEdge,2),:);
    for i = 1:size(V1,1)
        plot([V1(i,1) V2(i,1)],[V1(i,2) V2(i,2)],'k');
    end
	hold off          

    axis equal

    axis([-L L -L L])
    axis off
    set(gcf,'Color','white');
end