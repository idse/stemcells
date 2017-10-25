function instruct = makeInStruct(results, idx, debugplots, legendstr, dataDir, label)
    
    instruct = struct(  'A',[],'errA',[],'k',[],'errk',[],...
                    'Nmovies',0, 'Ncells',0,...
                    'ncrinit',[],'ncrfin',[],...
                    'nr',[],'cr',[],...
                    'Ninit',[],'Nfin',[],...
                    'Cinit',[],'Cfin',[],...
                    'bg',[],'x',[]);

    if debugplots
        h = [];
        h(1) = figure; hold on
        h(2) = figure; hold on
        h(3) = figure; hold on
        h(4) = figure; hold on
    end
    
    for i = 1:size(idx,1)

        result = results{idx(i,1)}{idx(i,2)};
        goodidx = logical(result.good);
        ncrstruct = getNCR(result);

        instruct.Nmovies = instruct.Nmovies + 1;
        instruct.Ncells = instruct.Ncells + sum(goodidx);

        instruct.A = cat(1, instruct.A, result.A(goodidx,1));
        instruct.errA = cat(1, instruct.errA, result.A(goodidx,1)-result.A(goodidx,2));
        instruct.k = cat(1, instruct.k, result.k(goodidx,1));
        instruct.errk = cat(1, instruct.errk, result.k(goodidx,1)-result.k(goodidx,2));
        
        % correct amplitude for cytoplasmic bleaching
        x = result.cytobleachFactor;
        if numel(x) == numel(goodidx)
            x = result.cytobleachFactor(goodidx);%(goodidx,1);
        else
            warning(['seg for ' result.description ' not up to date']);
        end
        instruct.Acorr = cat(1, instruct.A, result.A(goodidx,1)./x);
        
        instruct.ncrinit = cat(1, instruct.ncrinit, ncrstruct.ncrinit);
        instruct.ncrfin = cat(1, instruct.ncrfin, ncrstruct.ncrfin);
        
        instruct.nr = cat(1, instruct.nr, ncrstruct.nr);
        instruct.cr = cat(1, instruct.cr, ncrstruct.cr);
        
        instruct.Ninit = cat(1, instruct.Ninit, ncrstruct.Ninit);
        instruct.Nfin = cat(1, instruct.Nfin, ncrstruct.Nfin);
        instruct.Cinit = cat(1, instruct.Cinit, ncrstruct.Cinit);
        instruct.Cfin = cat(1, instruct.Cfin, ncrstruct.Cfin);
        instruct.bg = cat(1, instruct.bg, ncrstruct.bg);
        
        instruct.x = cat(1, instruct.x, ncrstruct.x);
        
        if debugplots
            xr = [ncrstruct.ncrinit]';
            figure(h(1))
            scatter(0*xr + i, xr, 500,'.');
            
            xr = [ncrstruct.ncrfin]';
            figure(h(2))
            scatter(0*xr + i, xr, 500,'.');

            figure(h(3))
            scatter(0*result.A(goodidx,1) + i, result.A(goodidx,1),500,'.');

            figure(h(4))
            scatter(0*result.k(goodidx,1) + i, result.k(goodidx,1),500,'.');
        end
    end
    if debugplots
        varnames = {'rin','rfin','A','k'};
        ylimvals = {[0 1.2],[0 1.2],[0 0.6],[0 0.025]};
        for i = 1:numel(h)
            figure(h(i))
            title(varnames(i));
            legend(legendstr,'Interpreter','none')
            xticks(1:size(idx,1));
            xlim([0.5 size(idx,1)+0.5]);
            ylim(ylimvals{i});
            hold off
            saveas(h(i), fullfile(dataDir,['debug_' varnames{i} label '.png']));
            saveas(h(i), fullfile(dataDir,['debug_' varnames{i} label '.fig']));
        end
    end
end