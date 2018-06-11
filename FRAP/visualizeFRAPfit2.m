function visualizeFRAPfit2(results,labels,onlygood)

    if ~exist('labels','var')
        labels = [];
    end
    
    if ~exist('onlygood','var')
        onlygood = false;
    end

    if strcmp(results.fitType, results.bleachType)
            func = @(p1,p2,p3,x) p1*(1-exp(-x*p2)) + p3;
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

    Nfrapped = size(tracesnorm,1);
    if onlygood && isfield(results,'good')
        shapesIdx = find(results.good);
    else
        shapesIdx = 1:Nfrapped;
    end
    
    t = results.tres*(0:size(tracesnorm,2)-1);
    tmax = results.tmax;
    tlim = round(max(tmax)*results.tres);
    frapframe = results.frapframe;

    lw = 2;
    colors = lines(Nfrapped);

    hold on 

    for shapeIdx = shapesIdx
        plot(t, tracesnorm(shapeIdx,:),...
                        'Color',colors(shapeIdx,:),'LineWidth',lw)
    end

    legendstr = {};
    for shapeIdx = shapesIdx
        A = results.A(shapeIdx,1);
        k = results.k(shapeIdx,1);
        B = results.B(shapeIdx,1);
            
        kerr = results.k(shapeIdx,3)-results.k(shapeIdx,1);
        Aerr = results.A(shapeIdx,3)-results.A(shapeIdx,1);
        Berr = results.B(shapeIdx,3) - results.B(shapeIdx,1);
        
        fitcurve = func(A,k,B,t(frapframe:tmax(shapeIdx))-t(frapframe));
        tau=1/(60*k);
        tauerr = kerr/(60*k^2);
        
        legendstr{shapeIdx} = ['A = ' num2str(A,2)...% '(' num2str(Aerr,1) ')'...
                                ';\tau=' num2str(tau,2) ';k=' num2str(k,2)];% '(' num2str(tauerr,1) ') min'];
        legendstr{shapeIdx} = [legendstr{shapeIdx} ', B = ' num2str(B,2)];
        
        plot(t(frapframe:tmax(shapeIdx)),fitcurve,...
                        'Color', colors(shapeIdx,:),'LineWidth',lw)
    end
    if ~isempty(labels)
        for shapeIdx = shapesIdx
            legendstr{shapeIdx} = [labels{shapeIdx} ': ' legendstr{shapeIdx}];
        end
    end

    hold off
    xlim([0 t(end)]);
    ylim([0 1.5*max(results.A(shapesIdx,1)+results.B(shapesIdx,1))]); 

    legend(legendstr(shapesIdx), 'Location','SouthEast');

    fs = 20;
    xlabel('time after bleach (sec)', 'FontSize',fs, 'FontWeight','Bold');
    ylabel('recovered fraction N^l/N_0', 'FontSize',fs, 'FontWeight','Bold');
    set(gcf,'color','w');
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
end
