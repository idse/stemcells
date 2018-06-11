function stress = computeDilationStress(g)

    nCells = numel(g.cells);
    stress = zeros([nCells 1]);

    for ci = 1:nCells

        vidx = g.bonds(g.cells{ci},1);
        vidx2 = g.bonds(g.cells{ci},2);
        verts = g.verts(vidx,:);

        edgevec = g.verts(vidx2,:) - g.verts(vidx,:);

        bonds = g.bonds;
        edgeBonds = bonds(bonds(:,4) == 0,:);
        dualEdgeBonds = edgeBonds(:,[2 1 4 3]);
        bonds = cat(1,bonds,dualEdgeBonds);

        outbondsPre = bonds(bonds(:,3) ~= ci & bonds(:,4) ~= ci,:);
        outbonds = zeros([size(g.verts,1) 4]);
        outbonds(outbondsPre(:,1),:) = outbondsPre;
        outbonds = outbonds(vidx,:);
        bidx = outbonds(:,1) > 0;
        outbonds = outbonds(bidx,:);

        outvec = g.verts(outbonds(:,2),:) - g.verts(outbonds(:,1),:);% 
        % hold on
        % 
        % scatter(g.verts(vidx(2:5),1),g.verts(vidx(2:5),2),'or');
        % 
        % quiver(verts(:,1),verts(:,2),outvec(:,1),outvec(:,2),0,'g','LineWidth',2);
        % quiver(verts(:,1),verts(:,2),edgevec(:,1),edgevec(:,2),0,'r','LineWidth',2);
        % 
        % hold off
        % axis off

        % normalize vectors
        edgevec = bsxfun(@rdivide, edgevec, sqrt(sum(edgevec.^2,2)));
        outvec = bsxfun(@rdivide, outvec, sqrt(sum(outvec.^2,2)));
        edgevecshift = circshift(edgevec,[1 0]);

        % sines of angles
        sinip = cross(outvec, edgevec(bidx,:));
        sinim = cross(outvec, edgevecshift(bidx,:));
        sinip = sinip(:,3);
        sinim = sinim(:,3);

        factor = sinim./circshift(sinip, [-1 0]);

        dilation = bsxfun(@times, outvec, factor);

        sigma = sum(sum(g.cellForce{ci}(bidx,:).*dilation(:,1:2), 2))/g.A(ci);
        if isnan(sigma) 
            sigma = 0;
        end
        stress(ci) = sigma;

%         isogonal = 0*verts;
%         isogonal(bidx, :) = 0.1*dilation;
%         isogonal = isogonal + verts;
    %     fv = struct();
    %     fv.Vertices = isogonal;
    %     fv.Faces = 1:6;
    %     patch(fv, 'FaceColor','none')
    end
end