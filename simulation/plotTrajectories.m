function plotTrajectories(X,t,params,b)
x = -2:0.1:2;
fs = 15;

plot(X',t)
hold off
xlim([x(1) x(end)])
xlabel('space');
ylabel('time');
legend('landscape', 'Location', 'SouthWest')
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')

if exist('params','var') && exist('b','var')
    hold on;
    V = @(x) polyval(params,x);
    plot(x,V(x)+b.*x, '-k','LineWidth',2);
    hold off
end