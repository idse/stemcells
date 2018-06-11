function plotFRAPcurves(dataDir, oibfile, results)

    tracesNuc = results.tracesNuc;
    tracesNucNorm = results.tracesNucNorm;

    tcut = size(tracesNucNorm,2);
    Nfrapped = size(tracesNucNorm,1);
    tres = results.tres;
    colors = hsv(Nfrapped);

    % filename
    [~,barefname,~] = fileparts(oibfile);
    barefname = strrep(barefname,'.','dot');

    % FRAP curves raw nuclear
    figure,
    hold on
    t = (1:tcut)*tres;
    legendstrs = {};
    for i = 1:Nfrapped
        legendstrs = [legendstrs num2str(i)];
        plot(t' ,tracesNuc(i,:)','Color',colors(i,:),'LineWidth',1.5);
    end
    hold off
    xlabel('time (sec)');
    ylabel('intensity')
    xlim([0 tcut*tres]);
    legend(legendstrs);
    saveas(gcf,fullfile(dataDir, [barefname '_FRAPcurvesNucRaw']));
    saveas(gcf,fullfile(dataDir, [barefname '_FRAPcurvesNucRaw'  '.png']));
    close;
    
    % FRAP curves raw cytoplasmic
    if isfield(results,'tracesCyt')
        figure,
        hold on
        for i = 1:Nfrapped
            plot(t' ,results.tracesCyt(i,:)','Color',colors(i,:),'LineWidth',1.5);
        end
        hold off
        xlabel('time (sec)');
        ylabel('intensity')
        xlim([0 tcut*tres]);
        legend(legendstrs);
        saveas(gcf,fullfile(dataDir, [barefname '_FRAPcurvesCytRaw']));
        saveas(gcf,fullfile(dataDir, [barefname '_FRAPcurvesCytRaw'  '.png']));
        close;
    end
    
    % FRAP curves normalized
    clf
    hold on
    for i = 1:Nfrapped
        plot(t' ,tracesNucNorm(i,:)','Color',colors(i,:),'LineWidth',1.5);
    end
    hold off
    xlabel('time (sec)');
    ylabel('normalized intensity')
    xlim([0 tcut*tres]);
    legend(legendstrs);
    saveas(gcf,fullfile(dataDir, [barefname '_FRAPcurvesNucNorm']));
    saveas(gcf,fullfile(dataDir, [barefname '_FRAPcurvesNucNorm.png']));
    
    if isfield(results,'tracesCyt')
        clf
        hold on
        for i = 1:Nfrapped
            plot(t' ,results.tracesCytNorm(i,:)','Color',colors(i,:),'LineWidth',1.5);
        end
        hold off
        xlabel('time (sec)');
        ylabel('normalized intensity')
        xlim([0 tcut*tres]);
        legend(legendstrs);
        saveas(gcf,fullfile(dataDir, [barefname '_FRAPcurvesCytNorm']));
        saveas(gcf,fullfile(dataDir, [barefname '_FRAPcurvesCytNorm.png']));
    end
end