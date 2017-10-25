function outstruct = getNCR(result)
    % outstruct = getNCR(result)

    out = [];
    outstruct = struct('x',[],'ncrinit',[],'ncrfin',[],...
                        'nr',[],'cr',[],...
                        'Ninit',[],'Nfin',[],'Cinit',[],'Cfin',[],'bg',[]);
    
	k = 1;
    for j = find(result.good)
        
        nucBeforeBleach = mean(result.tracesNuc(j,1:result.frapframe-1));
        %nucAfterBleach = mean(result.tracesNuc(j,result.frapframe:result.frapframe+1));
        nucAfterRecovery = mean(result.tracesNuc(j,end-10:end));
        tm = result.tmax(j);

        if isfield(result, 'tracesCyt')
            cytBeforeBleach = mean(result.tracesCyt(j,1:result.frapframe-1));
            cytAfterBleach = mean(result.tracesCyt(j,result.frapframe:result.frapframe+1));
            cytAfterRecovery = mean(result.tracesCyt(j,tm-10:tm));
            %r = (result.tracesNuc(j,:) - bg)./(result.tracesCyt(j,:) - bg);
        
        % below is just to deal with old formats without first reanalyzing
        % everything
        elseif ~iscell(result.cytstart)
            cytBeforeBleach = mean(result.cytstart(j,1:result.frapframe-1));
            cytAfterBleach = mean(result.cytstart(j,result.frapframe:result.frapframe+1));
            cytAfterRecovery = mean(result.cytend(j,tm-10:tm));
        else
            cytBeforeBleach = result.cytstart{j};
            cytAfterBleach = result.cytstart{j};
            cytAfterRecovery = result.cytend{j};
        end

        bg = 160;%min(result.tracesNuc(:,result.frapframe))%nucAfterBleach;
        x = (cytAfterBleach - bg)/(cytBeforeBleach - bg);
        ncrinit = (nucBeforeBleach - bg)/(cytBeforeBleach - bg);
        ncrfin = (nucAfterRecovery - bg)/(cytAfterRecovery - bg); 
        
        nr = (nucAfterRecovery - bg)/(nucBeforeBleach - bg);
        cr = (cytAfterRecovery - bg)/(cytBeforeBleach - bg); 
        
        outstruct(k) = struct('x',x,'ncrinit',ncrinit,'ncrfin',ncrfin,...
                                'nr',nr,'cr',cr,...
                        'Ninit',nucBeforeBleach,'Nfin',nucAfterRecovery,...
                        'Cinit',cytBeforeBleach,'Cfin',cytAfterRecovery,...
                        'bg', bg);
        k = k+1;
    end
end