function Ttot = meanTension(g)

    % add up tension of bond and antibond
    % good for the perimeter term

    bondstmp = g.bonds;
    Ttot = g.T;
    
    for i = 1:size(bondstmp,1)

         antibond = bondstmp(:,1) == bondstmp(i,2) & bondstmp(:,2) == bondstmp(i,1);

         if any(antibond)
            Ttot(i) = (g.T(i) + g.T(antibond))/2;
            Ttot(antibond) = (g.T(i) + g.T(antibond))/2;
            bondstmp(antibond,1) = 0;
         end

         bondstmp(i,1) = 0;
    end
end
    