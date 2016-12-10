function visualizeTimetrace(positions, meta, wellsWanted, saveDir)

    if nargin <= 3
        saveResult = false;
    else
        saveResult = true;
    end

    nWells = meta.nWells;
    posPerCondition = meta.posPerCondition;
    minNCells = 10; % minimal number of cells

    % plot settings
    fgc = 'k';
    bgc = 'w';
    graphbgc = 1*[1 1 1]; 
    graphfgc = 'r';
   
    s = strsplit(meta.timeInterval,' ');
    dt = str2double(s{1});
    unit = s{2};
    t = ((1:tmax) - treatmentTime)*dt;
    xlims = [t(1), t(end)+50];

    colors = hsv(posPerCondition); %0.5*[1 0 0]
    baseline = zeros([1 nWells]); % store baseline avg of each well

    for wellnr = wellsWanted

        conditionPositions = posPerCondition*(wellnr-1)+1:posPerCondition*wellnr;

        ttraceCat = cat(1,positions.timeTraces);
        ttraceCat = ttraceCat(conditionPositions);

        % weigh by number of cells, tends to make little difference
        W = cat(1,positions.ncells); 
        W = W(conditionPositions,:)';
        %W = ones(size(W)); % don't weigh
        W = bsxfun(@rdivide, posPerCondition*W, sum(W,2));

        nucTrace = cat(2,ttraceCat.nucLevelAvg);
        cytTrace = cat(2,ttraceCat.cytLevelAvg);
        bgTrace = cat(2,ttraceCat.background);

        nucMean = nanmean(nucTrace(1:tmax,:).*W,2);
        cytMean = nanmean(cytTrace(1:tmax,:).*W,2);
        bgMean = nanmean(bgTrace(1:tmax,:).*W,2);

        ratioMean = (nucMean - bgMean)./(cytMean - bgMean);
        %meanRatio = nanmean(ratio,1); % makes little difference
        baseline(wellnr) = mean(ratioMean(t < treatmentTime));

        clf
        hold on
        ratio = zeros([numel(conditionPositions) tmax]);
        for i = 1:numel(conditionPositions)

            pi = conditionPositions(i);

            nucTrace = positions(pi).timeTraces.nucLevelAvg;
            bgTrace = positions(pi).timeTraces.background;
            cytTrace = positions(pi).timeTraces.cytLevelAvg;

            R = (nucTrace - bgTrace)./(cytTrace - bgTrace);
            ratio(pi,:) = R(1:tmax)';
            ratio(pi, positions(pi).ncells < minNCells) = NaN;
            plot(t,ratio(pi,:),'Color', colors(i,:))
        end

        plot(t, ratioMean, graphfgc,'LineWidth',2)
        %plot(t, meanRatio, 'g','LineWidth',2)

        legend({'1','2','3','4','mean'});

        fs = 24;
        xlabel(['time (' unit ')'], 'FontSize',fs,'FontWeight','Bold','Color',fgc)
        ylabel('nuclear : cytoplasmic Smad4', 'FontSize',fs,'FontWeight','Bold','Color',fgc);

        xlim(xlims);
        ylim(ylims);
        
        set(gcf,'color',bgc);
        set(gca, 'LineWidth', 2);
        set(gca,'FontSize', fs)
        set(gca,'FontWeight', 'bold')
        set(gca,'XColor',fgc);
        set(gca,'YColor',fgc);
        set(gca,'Color',graphbgc);

        if saveResult
            export_fig(fullfile(dataDir, ['timeTrace_well' num2str(wellnr) '.png']),'-native -m2');
        end

        hold off
    end
end