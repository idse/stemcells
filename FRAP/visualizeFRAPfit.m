function visualizeFRAPfit(results)

    func = @(p1,p2,x) p1*(1-exp(-x*p2));

    if strcmp(results.bleachType,'cytoplasmic')
        traces = results.tracesCyt;
        tracesnorm = results.tracesCytNorm;
    else
        traces = results.tracesNuc;
    end
    
    t = results.tres*(0:size(traces,2)-1);
    tmax = results.tmax;
    tlim = round(max(tmax)*results.tres);
    frapframe = results.frapframe;
    
    Nfrapped = size(traces,1);
    
    lw = 2;
    colors = lines(Nfrapped);
    
    clf
    hold on 
    
    for shapeIdx = 1:Nfrapped
        plot(t, tracesnorm(shapeIdx,:),...
                        'Color',colors(shapeIdx,:),'LineWidth',lw)
    end

    legendstr = {};
    for shapeIdx = 1:Nfrapped
        A = results.A(shapeIdx,1);
        k = results.k(shapeIdx,1);
        kerr = results.k(shapeIdx,3)-results.k(shapeIdx,1);
        Aerr = results.A(shapeIdx,3)-results.A(shapeIdx,1);
        tau=1/(60*k);
        tauerr = kerr/(60*k^2);
        plot(t(frapframe:tmax(shapeIdx)),func(A,k,t(frapframe:tmax(shapeIdx))),...
                        'Color', colors(shapeIdx,:),'LineWidth',lw)
        legendstr{shapeIdx} = ['A = ' num2str(A,2) '(' num2str(Aerr,1)...
                            '), \tau=' num2str(tau,2) '(' num2str(tauerr,1) ') min'];
    end

    hold off
    xlim([0 t(end)]);
    ylim([0 1.2*max(results.A(:,1))]); 
    legend(legendstr, 'Location','SouthEast');

    fs = 20;
    xlabel('time after bleach (sec)', 'FontSize',fs, 'FontWeight','Bold');
    ylabel('recovered fraction N^l/N_0', 'FontSize',fs, 'FontWeight','Bold');
    set(gcf,'color','w');
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    
    % Aryeh's version of the plot:
    % subplot(1,2,shapeIdx); plot(tdata,fdata,'r.'); hold on; plot(outfit{shapeIdx},'k','predfunc');
end
