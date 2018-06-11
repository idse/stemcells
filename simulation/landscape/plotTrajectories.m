function plotTrajectories(X,t,params,b)
% Plot a set of trajectories and (optionally) a landscape
% X are the trajectories, rows are cells and columns timepoints
% params defines the polynomial for the landscape and b the biasing force.
% if both params and b are supplied, the landscape will be plotted together
% with the trajectories. Landscape is auto-scaled between tmin and tmax.
% see also runTrajectories

x = -2:0.1:2;
fs = 36;

plot(X',t)
hold off
xlim([x(1) x(end)])
xlabel('space');
ylabel('time');
set(gcf,'color','w');
set(gca, 'LineWidth', 3);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')

if exist('params','var') && exist('b','var')
    hold on;
    %compute landscape
    lscape =  polyval(params,x)+b.*x;
    
    %rescale the landscape for y values of the plot
    tmin = min(t); tmax = max(t);
    lmin = min(lscape); 
    
    lscape = lscape+(tmin-lmin);
    lmax = max(lscape);
    lscape = lscape*tmax/lmax;
    
    
    plot(x,lscape, '-k','LineWidth',2);
    hold off
end