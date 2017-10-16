function out = getNCR(result)
    % out = getNCR(result)
    % out ~ [x, rinit, rfin]

    out = [];
    
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
        rinit = (nucBeforeBleach - bg)/(cytBeforeBleach - bg);
        rfin = (nucAfterRecovery - bg)/(cytAfterRecovery - bg); 
        
        out = cat(1, out, [x, rinit, rfin]);
    end
end