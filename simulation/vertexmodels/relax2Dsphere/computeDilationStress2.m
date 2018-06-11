function stress = computeDilationStress2(g)
    % using normal dilation
    
    nCells = numel(g.cells);
    stress = zeros([nCells 1]);

    for ci = 1:nCells

        vidx = g.bonds(g.cells{ci},1);
        verts = g.verts(vidx,:);
        vertsCM = bsxfun(@minus, verts, mean(verts));
        dilation = bsxfun(@rdivide, vertsCM, sqrt(sum(vertsCM.^2,2)));
        
        sigma = sum(sum(g.cellForce{ci}.*dilation(:,1:2), 2))/g.A(ci);
        if isnan(sigma) 
            sigma = 0;
        end
        stress(ci) = sigma;
    end
end