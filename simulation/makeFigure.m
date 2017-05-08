params=  [0.25 0 -1 0 0];
b =1;
[X, t, Rfrac]=runTrajectories(params,b);
plotTrajectories(X,t,params,b);
ylim([min(t) max(t)]);
export_fig('~/Dropbox/Grants/Simons2016/Figures/1a.eps','-transparent');

params=  [0.25 0 -1 0 0];
b =0;
[X, t, Rfrac]=runTrajectories(params,b);
plotTrajectories(X,t,params,b);
ylim([min(t) max(t)]);

export_fig('~/Dropbox/Grants/Simons2016/Figures/1b.eps','-transparent');

params=  [0.25 0 -1 0 0];
b =-1;
[X, t, Rfrac]=runTrajectories(params,b);
plotTrajectories(X,t,params,b);
ylim([min(t) max(t)]);

export_fig('~/Dropbox/Grants/Simons2016/Figures/1c.eps','-transparent');