function gnew = T1transition(g, bi)
    % do a T1 transition
    % 
    % gnew = T1transition(g, bi)
    %
    % bi:   bond index

    v1i = g.bonds(bi,1);
    v2i = g.bonds(bi,2);
    c1i = g.bonds(bi,3);
    c2i = g.bonds(bi,4);
    
    if ~(c1i && c2i)
        warning('no T1s are allowed on the boundary, returning');
        gnew = g;
        return
    end
    
    % antibond index 
    abi = find(g.bonds(:,3) == g.bonds(bi,4) & g.bonds(:,4) == g.bonds(bi,3));

    % BELOW WILL NOT WORK IF THERE WERE 4-FOLDS
    
    % other cells involved
    bintov1 = g.bonds(:,2) == v1i & g.bonds(:,3) == c1i;
    c3i = g.bonds(bintov1, 4);
    bintov2 = g.bonds(:,1) == v2i & g.bonds(:,3) == c1i;
    c4i = g.bonds(bintov2, 4);
    
    if isempty(c3i), c3i = 0; end    
    if isempty(c4i), c4i = 0; end
    
    disp(['T1 on cells ' num2str([c1i c2i c3i c4i])]);
    
    % bond, previous and next within the cell
	cbi = g.cells{c1i} == bi;
    bin = g.cells{c1i}(circshift(cbi, [0 1]));
    bip = g.cells{c1i}(circshift(cbi, [0 -1]));
    
    % antibond, previous and next within the cell
    cabi = g.cells{c2i} == abi;
    abin = g.cells{c2i}(circshift(cabi, [0 1]));
    abip = g.cells{c2i}(circshift(cabi, [0 -1]));
  
    % delete bonds and antibond from c1 and c2
    g.cells{c1i}(g.cells{c1i} == bi) = [];
    g.bonds(bin,1) = g.bonds(bip,2);

    g.cells{c2i}(g.cells{c2i} == abi) = [];
    g.bonds(abin,1) = g.bonds(abip,2);
    
    % then add them to c3 and c4, if they are not the outside cell
    if c3i
        
        b3 = find(g.bonds(:,2) == v1i & g.bonds(:,3) == c3i);
        cb3 = find(g.cells{c3i} == b3);

        g.bonds(g.cells{c3i}(cb3),2) = v2i;
        
        if cb3 < numel(g.cells{c3i})
            g.bonds(g.cells{c3i}(cb3+1),1) = v1i;
            g.cells{c3i} = [g.cells{c3i}(1:cb3) abi g.cells{c3i}(cb3+1:end)];
        else
            g.bonds(g.cells{c3i}(1),1) = v1i;
            g.cells{c3i} = [g.cells{c3i}(1:cb3) abi];
        end
            
        g.bonds(abi,3:4) = [c3i c4i];
    else
        g.bonds(abi,:) = [0 0 0 0];
    end
    
    % and c4
    if c4i
        
        b4 = find(g.bonds(:,2) == v2i & g.bonds(:,3) == c4i);
        cb4 = find(g.cells{c4i} == b4);

        g.bonds(g.cells{c4i}(cb4),2) = v1i;
        g.bonds(g.cells{c4i}(cb4+1),1) = v2i;
        g.cells{c4i} = [g.cells{c4i}(1:cb4) bi g.cells{c4i}(cb4+1:end)];

        g.bonds(bi,3:4) = [c4i c3i];
    else
        g.bonds(bi,:) = [0 0 0 0];
    end
    
    gnew = g;
end