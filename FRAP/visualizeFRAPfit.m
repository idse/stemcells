function visualizeFRAPfit(results,labels)

    if ~exist('labels','var')
        labels = [];
    end

    if strcmp(results.fitType, results.bleachType)
            func = @(p1,p2,x) p1*(1-exp(-x*p2));
            decay = false;
    else 
            func = @(p1,p2,p3,x) p1*exp(-x*p2) + p3;
            decay = true;
    end
    
    if strcmp(results.fitType,'cytoplasmic')
        tracesnorm = results.tracesCytNorm;
    else
        tracesnorm = results.tracesNucNorm;
    end

    t = results.tres*(0:size(tracesnorm,2)-1);
    tmax = results.tmax;
    tlim = round(max(tmax)*results.tres);
    frapframe = results.frapframe;
    
    Nfrapped = size(tracesnorm,1);
    
    lw = 2;
    colors = lines(Nfrapped);

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
        
        legendstr{shapeIdx} = ['A = ' num2str(A,2)...% '(' num2str(Aerr,1) ')'...
                                ', \tau=' num2str(tau,2) ' min'];% '(' num2str(tauerr,1) ') min'];
        
        if ~decay
            fitcurve = func(A,k,t(frapframe:tmax(shapeIdx)));
        else
            B = results.B(shapeIdx,1);
            Berr = results.B(shapeIdx,3) - results.B(shapeIdx,1);
            fitcurve = func(A,k,B,t(frapframe:tmax(shapeIdx)));
            legendstr{shapeIdx} = [legendstr{shapeIdx} ', B = ' num2str(B,2)];
        end
        plot(t(frapframe:tmax(shapeIdx)),fitcurve,...
                        'Color', colors(shapeIdx,:),'LineWidth',lw)
    end
    if ~isempty(labels)
        for shapeIdx = 1:Nfrapped
            legendstr{shapeIdx} = [labels{shapeIdx} ': ' legendstr{shapeIdx}];
        end
    end

    hold off
    xlim([0 t(end)]);
    if decay
        ylim([0 1.2*max(results.A(:,1)+results.B(:,1))]); 
    else
        ylim([0 1.2*max(results.A(:,1))]); 
    end
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
