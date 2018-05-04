function newg = divideCell(gcur, ci, divaxis)

    offset = 0;

    [faces,bonds,verts, nfi, nbi, ibi] = divideFace(gcur.cells, gcur.bonds,...
                                                    gcur.verts,...
                                                    ci, divaxis, offset);

    newg = gcur;
    newg.cells = faces;
    newg.bonds = bonds;
    newg.verts = verts;
    
%     % new preferred area is inherited from mother
%     newg.A0 = cat(1,newg.A0, gcur.A0(ci));

    % new preferred area coefficient is a quarter of the mother
    % in the E ~ A0/A model
    gcur.A0(ci) =  gcur.A0(ci)/(2*sqrt(2));
    % in the preferred area model (A-A0)^2, we have half the area, twice
    % the k
    %gcur.A0(ci) =  gcur.A0(ci)/2;
    newg.A0 = cat(1,gcur.A0, gcur.A0(ci));
    
    %gcur.kA0(ci) = 2*gcur.kA0(ci);
    %newg.kA0 = cat(1,newg.kA0, gcur.kA0(ci));
    
    % new preferred perimeter is also inherited from mother
    newg.p0 = cat(1,newg.p0, gcur.p0(ci));
    
    % as is the perimeter tension
    gcur.pT0(ci) = gcur.pT0(ci);
    newg.pT0 = cat(1,gcur.pT0, gcur.pT0(ci));
    
    % and the stress to prevent stress response in dividing cells
    newg.stress = cat(1,newg.stress, gcur.stress(ci));

    % new tension is inherited too:
    % intersected bonds will have the same tension on each side 
    % the partition itself will have the mean of the bonds it intersects
    Told0 = newg.T0(ibi);
    newg.T0 = cat(1, newg.T0, ones([6 1]));
    newg.T0(nbi(1:2)) = mean(Told0);
    newg.T0(nbi(3:4)) = Told0(1);
    newg.T0(nbi(5:6)) = Told0(2);
    
    % preferred length the same way
    l0old = newg.l0(ibi);
    newg.l0 = cat(1, newg.l0, ones([6 1]));
    newg.l0(nbi(1:2)) = mean(l0old);
    newg.l0(nbi(3:4)) = l0old(1);
    newg.l0(nbi(5:6)) = l0old(2);
    
    % inherit clonal identity;
    newg.clones = cat(1, newg.clones, gcur.clones(ci));
end