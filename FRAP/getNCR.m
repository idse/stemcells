function outstruct = getNCR(result, bleachcorrect)
    % outstruct = getNCR(result)

    outstruct = struct('x',[],'ncrinit',[],'ncrfin',[],...
                        'nr',[],'cr',[],...
                        'Ninit',[],'Nbleach',[],'Nfin',[],...
                        'Cinit',[],'Cbleach',[],'Cfin',[],'bg',[]);
	k = 1;

    if ~exist('bleachcorrect','var')
        bleachcorrect = false;
    end

    if bleachcorrect
        tracesNuc = result.tracesNucBC;
        tracesCyt = result.tracesCytBC;
    else
        tracesNuc = result.tracesNuc;
        tracesCyt = result.tracesCyt;
    end
    
    for j = find(result.good)
        
        tm = result.tmax(j);
        %tm = size(result.tracesNuc,2)*ones([size(result.tracesNuc,1) 1]);

        nucBeforeBleach = nanmean(tracesNuc(j,1:result.frapframe-1));
        nucAfterBleach = nanmean(tracesNuc(j,result.frapframe:result.frapframe+1));
        %nucAfterRecovery = mean(result.tracesNuc(j,tm-10:tm));
        nucAfterRecovery = nanmean(tracesNuc(j,tm-10:tm));

        if isfield(result, 'tracesCyt')
            cytBeforeBleach = nanmean(tracesCyt(j,1:result.frapframe-1));
            % for cyt : slightly later postbleach window to get rid of
            % diffusive recovery
            %cytAfterBleach = nanmean(tracesCyt(j,result.frapframe+2:result.frapframe+5));
            [~,i] = max(tracesCyt(j,result.frapframe+2:result.frapframe+20));
            cytAfterBleach = nanmean(tracesCyt(j,result.frapframe+i:result.frapframe+i+2));
            
            cytAfterRecovery = nanmean(tracesCyt(j,tm-10:tm));
            %r = (result.tracesNuc(j,:) - bg)./(result.tracesCyt(j,:) - bg);
        
        % below is just to deal with old formats without first reanalyzing
        % everything
        elseif ~iscell(result.cytstart)
            cytBeforeBleach = nanmean(result.cytstart(j,1:result.frapframe-1));
            cytAfterBleach = nanmean(result.cytstart(j,result.frapframe+1:result.frapframe+3));
            cytAfterRecovery = nanmean(result.cytend(j,tm-10:tm));
        else
            cytBeforeBleach = result.cytstart{j};
            cytAfterBleach = result.cytstart{j};
            cytAfterRecovery = result.cytend{j};
        end

        %bg = 120;%min(result.tracesNuc(:,result.frapframe))%nucAfterBleach;
        %bg = min(min(result.tracesNuc(:,1:10)));
        %bg = min(bg, result.bgempty);
        bg = result.bg;
        
        x = (cytAfterBleach - bg)/(cytBeforeBleach - bg);
        ncrinit = (nucBeforeBleach - bg)/(cytBeforeBleach - bg);
        ncrfin = (nucAfterRecovery - bg)/(cytAfterRecovery - bg); 

        nr = (nucAfterRecovery - bg)/(nucBeforeBleach - bg);
        cr = (cytAfterRecovery - bg)/(cytBeforeBleach - bg); 
        
        outstruct(k) = struct('x',x,'ncrinit',ncrinit,'ncrfin',ncrfin,...
                                'nr',nr,'cr',cr,...
                        'Ninit',nucBeforeBleach,'Nbleach',nucAfterBleach,'Nfin',nucAfterRecovery,...
                        'Cinit',cytBeforeBleach,'Cbleach',cytAfterBleach','Cfin',cytAfterRecovery,...
                        'bg', bg);
        k = k+1;
    end
end