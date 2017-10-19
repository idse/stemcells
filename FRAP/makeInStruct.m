function instruct = makeInStruct(results, idx, legendstr, dataDir, label)
    
    instruct = struct(  'A',[],'errA',[],'k',[],'errk',[],...
                    'Nmovies',0, 'Ncells',0,...
                    'rinit',[],'rfin',[]);
                
    h = [];
    h(1) = figure; hold on
    h(2) = figure; hold on
    h(3) = figure; hold on
    h(4) = figure; hold on
    
    for i = 1:size(idx,1)

        result = results{idx(i,1)}{idx(i,2)};
        goodidx = logical(result.good);
        xr = getNCR(result);

        instruct.Nmovies = instruct.Nmovies + 1;
        instruct.Ncells = instruct.Ncells + sum(goodidx);

        instruct.rinit = cat(1, instruct.rinit, xr(:,2));
        instruct.rfin = cat(1, instruct.rfin, xr(:,3));
        instruct.A = cat(1, instruct.A, result.A(goodidx,1));
        instruct.errA = cat(1, instruct.errA, result.A(goodidx,1)-result.A(goodidx,2));
        instruct.k = cat(1, instruct.k, result.k(goodidx,1));
        instruct.errk = cat(1, instruct.errk, result.k(goodidx,1)-result.k(goodidx,2));
                
        figure(h(1))
        scatter(0*xr(:,2) + i, xr(:,2),500,'.');

        figure(h(2))
        scatter(0*xr(:,3) + i, xr(:,3),500,'.');

        figure(h(3))
        scatter(0*result.A(goodidx,1) + i, result.A(goodidx,1),500,'.');
        
        figure(h(4))
        scatter(0*result.k(goodidx,1) + i, result.k(goodidx,1),500,'.');
    end
    
    varnames = {'rin','rfin','A','k'};
    ylimvals = {[0 1.2],[0 1.2],[0 0.6],[0 0.03]};
    for i = 1:numel(h)
        figure(h(i))
        title(varnames(i));
        legend(legendstr,'Interpreter','none')
        xticks(1:size(idx,1));
        xlim([0.5 size(idx,1)+0.5]);
        ylim(ylimvals{i});
        hold off
        saveas(h(i), fullfile(dataDir,['debug_' varnames{i} label '.png']));
    end
end