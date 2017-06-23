function plotFRAPcurves(dataDir, oibfile, results)

    tracesnorm = results.tracesnorm;
    tcut = size(tracesnorm,2);
    Nfrapped = size(tracesnorm,1);
    tres = results.tres;
    
    tres
    % filename
    [~,barefname,~] = fileparts(oibfile);
    barefname = strrep(barefname,'.','dot');

    % FRAP curves
    figure,
    t = repmat((1:tcut)*tres,[Nfrapped 1]);
    plot(t' ,tracesnorm');
    xlabel('time (sec)');
    ylabel('intensity')
    saveas(gcf,fullfile(dataDir, ['FRAPcurvesRaw_' barefname]));
    saveas(gcf,fullfile(dataDir, ['FRAPcurvesRaw_' barefname '.png']));

    plot(t' ,tracesnorm');
    xlabel('time (sec)');
    ylabel('normalized intensity')
    saveas(gcf,fullfile(dataDir, ['FRAPcurvesNorm_' barefname]));
    saveas(gcf,fullfile(dataDir, ['FRAPcurvesNorm_' barefname '.png']));
    close;

    % FRAP fit
    visualizeFRAPfit(results)
    xlim([0 1600]);
    saveas(gcf,fullfile(dataDir, ['FRAPfit_' barefname]));
    saveas(gcf,fullfile(dataDir, ['FRAPfit_' barefname '.png']));
    close;
end