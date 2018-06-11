function plotDistribution(stats, channelIndex, cumulative)

    if ~exist('cumulative','var')
        cumulative = false;
    end
    fs = 15;
    
    nall = [stats.histograms{:, channelIndex}];
    dist = bsxfun(@rdivide, nall, sum(nall,1));
    if cumulative
        dist = cumsum(dist);
    end

    [x,y] = histForBarlikePlot(stats.bins{channelIndex}, dist);

    figure, 

    plot(x,y,'LineWidth',2);
    
    title(stats.channelLabel(channelIndex), 'FontSize',fs, 'FontWeight','Bold');

    xlabel('Intensity', 'FontSize',fs, 'FontWeight','Bold');
    ylabel('Frequency', 'FontSize',fs, 'FontWeight','Bold');
    
    set(gcf,'color','w');
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')

    xlim(stats.lim{channelIndex});
    
    if cumulative
        legend(stats.conditions, 'Location','SouthEast');
        ylim([0 1.1]);
    else
        legend(stats.conditions, 'Location','NorthEast');
    end
end