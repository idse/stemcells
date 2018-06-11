function [faces,bonds,verts, nfi, nbi, ibi] = divideFace(faces,bonds,verts,fi,daxis,fioffset)
% fi = ci;
% fioffset = nCells;
% faces = gfinal.facesA;
% bonds = gfinal.bonds;
% verts = gfinal.verts;
%
% nbi : new bond indices
% the first two are the partition
% the second two are daughters of ibi(1)
% the third two are daughters of ibi(2)
%
% ibi : the original bonds intersected by division
% nfi : new face indices

    % indices of new bonds, 6 bonds are created and added at the end of the
    % array (returning indices for generality)
    nbi = (size(bonds,1) + 1):(size(bonds,1) + 6);

    % vertex positions of the cell relative to CM
    vertOfCell = verts(bonds(faces{fi},1),1:2);
    cellCM = mean(vertOfCell);
    relVert = vertOfCell;
    relVert(:,1) = relVert(:,1) - cellCM(1);
    relVert(:,2) = relVert(:,2) - cellCM(2);

    % assign vertex indices to one cell or the other
    
    % logical index of bonds whose 1st vertex will be in daughter 1
    c1idx = false([1 numel(faces{fi})]);
    for i = 1:numel(faces{fi})
        if dot(relVert(i,:),daxis) > 0 
            c1idx(i) = true;
        end
    end
    
    % shift so that the first bond of daughter1 is first
    % this ensures c1bi(end) is the last bond ccw later 
    shift = numel(faces{fi}) - find(c1idx - circshift(c1idx,[0 -1]) < 0);
    faces{fi} = circshift(faces{fi}, [0 shift]);
    c1idx = circshift(c1idx, [0 shift]);
    
    % indices of bonds whose first vertex is in cell 1 or cell 2
    % the last bond in the list is the one intersected by division 
    c1bi = faces{fi}(c1idx);
    c2bi = faces{fi}(~c1idx);
    
    % for output
    ibi = [c1bi(end) c2bi(end)];
    
%     clf
%     patch(verts(bonds(faces{fi},1),1),verts(bonds(faces{fi},1),2),'w')
%     hold on
%     scatter(verts(bonds(c1bi,1),1),verts(bonds(c1bi,1),2),'g')
%     scatter(verts(bonds(c2bi,1),1),verts(bonds(c2bi,1),2),'r')
%     hold off
%     axis equal

    % create two additional vertices 
    newVert1 = mean(verts(bonds(c1bi(end),1:2),:));
    newVert2 = mean(verts(bonds(c2bi(end),1:2),:));
    nv1i = size(verts,1) + 1;
    nv2i = size(verts,1) + 2;
    verts = cat(1, verts, newVert1, newVert2);

%     hold on
%     scatter(newVert1(1),newVert1(2),'b')
%     scatter(newVert2(1),newVert2(2),'b')
%     hold off
    
    % new cell 1 will take place of mother in the table
    % index of the new cell 2
    nfi = numel(faces) + 1;

    % create new bonds:
    % 1 bond connecting both new vertices, partitioning the cells
    partition = [nv1i nv2i fi nfi];
    antipart = [nv2i nv1i nfi fi];
    pi = size(bonds,1) + 1;
    api = size(bonds,1) + 2;
    bonds = cat(1, bonds, partition, antipart);

    % 4 bonds replacing the old bonds that are intersected by the partition
    newBond1 = bonds(c2bi(end),:);
    newBond2 = bonds(c1bi(end),:);
    
    % the new bonds belong to a new cell, we update the cell index of all 
    % bonds in the new cell further down, when the cell is made
    newBond1(2) = nv2i;
    newBond2(1) = nv1i;
    
    bonds(c2bi(end),1) = nv2i;
    bonds(c1bi(end),2) = nv1i;
    nb1i = size(bonds,1) + 1;
    nb2i = size(bonds,1) + 2;
    bonds = cat(1, bonds, newBond1, newBond2);
    
    % two neighboring cells now have extra bond too
    % this is dividing up the anti bonds to the bonds that were just
    % partitioned
    
    % identify the neighbor
    nn1i = bonds(c2bi(end),4) - fioffset;
    nn2i = bonds(c1bi(end),4) - fioffset;

    if nn1i > 0
        
        % first update the antibond to the bond that already existed
        abidx1 = bonds(:,3) == bonds(c2bi(end),4) & bonds(:,4) == bonds(c2bi(end),3);
        bonds(abidx1,2) = bonds(c2bi(end),1);

        % now make the new antibond and add to table
        newAB1 = newBond1([2 1 4 3]);
        nab1i = size(bonds,1) + 1;
        bonds = cat(1, bonds, newAB1);
        
        % update the neighbor
        insertidx1 = find(faces{nn1i} == find(abidx1));
        nn1val = numel(faces{nn1i});
        faces{nn1i} = [faces{nn1i}(1:insertidx1) nab1i faces{nn1i}(insertidx1+1:nn1val)];
    end
    
    if nn2i > 0
        % first update the antibond to the bond that already existed
        abidx2 = bonds(:,3) == bonds(c1bi(end),4) & bonds(:,4) == bonds(c1bi(end),3);
        bonds(abidx2,1) = bonds(c1bi(end),2);

        % now make the new antibond and add to table
        newAB2 = newBond2([2 1 4 3]);
        nab2i = size(bonds,1) + 1;
        bonds = cat(1, bonds, newAB2);
        
        % update the neighbor
        insertidx2 = find(faces{nn2i} == find(abidx2));
        nn2val = numel(faces{nn2i});
        faces{nn2i} = [faces{nn2i}(1:insertidx2-1) nab2i faces{nn2i}(insertidx2:nn2val)];
    end
    
    % create the daughter cells
    faces{fi} = [c1bi pi c2bi(end)];
    faces{nfi} = [c2bi(1:end-1) nb1i api nb2i];
    
    % update bond and antibond face index for new face
    for i = 1:numel(faces{nfi})
        
        cbi = faces{nfi}(i);
        bonds(cbi,3) = nfi;
        cabi = bonds(:,2) == bonds(cbi,1) & bonds(:,1) == bonds(cbi,2);
        bonds(cabi,4) = nfi; 
    end
   
    %patch(verts(bonds(faces{fi},1),1),verts(bonds(faces{fi},1),2),'w')

    %nn = [nn1i nn2i];
    %newvert = [newVert1; newVert2];
end