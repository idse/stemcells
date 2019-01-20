function [kymograph, ybins] = makeKymograph(colonies, param)
    %
    % param : parameter structure with fields
    %   - xres    : xy resolution
    %   - tres    : time resolution
    %   - Irange  : intensity range
    %   - revbins : reverse bins boolean
    %   - rescale : rescale Irange from 0 to 1
    %   - label   : x-axis label
    %   - trange  : time range
    
    kymographs = {};

    Ninterp = 100;
    Ntime = numel(colonies(1).radialProfile);
    ci = 1;
    
    if ~isfield(param,'trange')
        param.trange = [1 Ntime];
    end
        
    for coli = 1:numel(colonies)

        if ~isempty(colonies(coli).radiusPixel)

        kymograph = NaN([Ntime Ninterp]);

        for ti = param.trange(1):param.trange(2)

            if ~isempty(colonies(coli).radialProfile(ti).NucAvg)

                bins = colonies(coli).radialProfile(ti).BinEdges(1:end-1)*param.xres';
                
                if param.revbins
                    bins = max(bins) - bins;
                end
                linbins = linspace(bins(1),bins(end),Ninterp);
                
                N = colonies(coli).radialProfile(ti).NucAvg(:,ci);
                C = colonies(coli).radialProfile(ti).CytAvg(:,ci);
                kymograph(ti,:) = interp1(bins, N./C, linbins);
            end
        end

        kymographs{coli} = kymograph;
        
        end
    end

	kymograph = nanmean(cat(3,kymographs{:}),3);
    
    if param.rescale
        kymograph(kymograph < param.Irange(1)) = param.Irange(1);
        kymograph(kymograph > param.Irange(2)) = param.Irange(2);
        kymograph = (kymograph - min(kymograph(:)))./(max(kymograph(:)) - min(kymograph(:)));
        param.Irange = [0 1];
    end

    tframe = 1:size(kymograph,1);
    thr = tframe*param.tres/60;
    ybins = linbins;

    figure,
    imagesc(kymograph','YData',ybins,'XData',thr, param.Irange); %[0.35 0.45] for N/(C+N)

    axis square
    fs = 28;
    
    if param.revbins
        ystr = 'edge';
    else
        ystr = 'center';
    end
    ylabel([ystr ' distance (um)    '], 'FontSize',fs, 'FontWeight','Bold');
    xlabel('time (hrs)', 'FontSize',fs, 'FontWeight','Bold');
    ylim([min(bins) max(bins)])
    
    set(gcf,'color','w');
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca,'Ydir','Normal')

    h = colorbar('Ticks',param.Irange);
    %xticks([0 5 10 15 20]);
    xlabel(h, param.label,'LineWidth', 2, 'FontSize', fs,'FontWeight', 'bold')
    %title('beta-catenin');
end