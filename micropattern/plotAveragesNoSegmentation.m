function plotAveragesNoSegmentation(datadir)

load(fullfile(datadir,'colonies.mat'));
load(fullfile(datadir,'metaData.mat'));
DAPIChannel = 1;

r = imfilter(colonies(1).radialProfile.BinEdges,[1 1]/2)*meta.xres;
r(1) = colonies(1).radialProfile.BinEdges(1)*meta.xres;
r = r(1:end-1);
colCat = cat(3,colonies(:).radialProfile);

nucAvgAll = mean(cat(3,colCat.NucAvg),3);

nucAvgAllNormalized = bsxfun(@rdivide, nucAvgAll, nucAvgAll(:,DAPIChannel));

plot(r, nucAvgAllNormalized(:,2:4),'.-','LineWidth',3)
legend(meta.channelLabel(2:4));
axis([min(r) max(r) 0 4]);

% how normalize cytoplasmic levels?%cytAvgAll = mean(cat(3,colCat.CytAvg),3);

%cytAvgAllNormalized = bsxfun(@rdivide, cytAvgAll, nucAvgAll(:,DAPIChannel));


% hold on
% %plot(r, cytAvgAllNormalized(:,2:4),'--')
% hold off