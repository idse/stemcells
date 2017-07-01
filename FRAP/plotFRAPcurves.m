function plotFRAPcurves(dataDir, oibfile, results)

    tracesNuc = results.tracesNuc;
    tracesNucNorm = results.tracesNucNorm;

    tcut = size(tracesNucNorm,2);
    Nfrapped = size(tracesNucNorm,1);
    tres = results.tres;
    colors = lines(Nfrapped);

    % filename
    [~,barefname,~] = fileparts(oibfile);
    barefname = strrep(barefname,'.','dot');

    % FRAP curves
    figure,
    hold on
    t = (1:tcut)*tres;
    for i = 1:Nfrapped
        plot(t' ,tracesNuc(i,:)','Color',colors(i,:),'LineWidth',1.5);
        if isfield(results,'tracesCyt')
            plot(t' ,results.tracesCyt(i,:)','Color',colors(i,:)*0.8,'LineWidth',1);
        end
    end
    hold off
    xlabel('time (sec)');
    ylabel('intensity')
    xlim([0 tcut*tres]);
    saveas(gcf,fullfile(dataDir, [barefname '_FRAPcurvesRaw']));
    saveas(gcf,fullfile(dataDir, [barefname '_FRAPcurvesRaw'  '.png']));
    close;
    
    clf
    hold on
    t = (1:tcut)*tres;
    for i = 1:Nfrapped
        plot(t' ,tracesNucNorm(i,:)','Color',colors(i,:),'LineWidth',1.5);
        if isfield(results,'tracesCyt')
            plot(t' ,results.tracesCytNorm(i,:)','Color',colors(i,:)*0.8,'LineWidth',1);
        end
    end
    hold off
    xlabel('time (sec)');
    ylabel('normalized intensity')
    xlim([0 tcut*tres]);
    saveas(gcf,fullfile(dataDir, [barefname '_FRAPcurvesNorm']));
    saveas(gcf,fullfile(dataDir, [barefname '_FRAPcurvesNorm.png']));

    % FRAP fit
    visualizeFRAPfit(results)
    xlim([0 tcut*tres]);
    saveas(gcf,fullfile(dataDir, [barefname '_FRAPfit' ]));
    saveas(gcf,fullfile(dataDir, [barefname '_FRAPfit.png']));
    close;
end