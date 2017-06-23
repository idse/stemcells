function [nucAvgAllNormalized, r] = plotAveragesNoSegmentation(datadir,colSize,DAPIChannel)

if ~exist('DAPIChannel','var')
    DAPIChannel = 1;
end

load(fullfile(datadir,'colonies.mat'));
load(fullfile(datadir,'metaData.mat'));

if ~exist('colSize','var')
    colSize = colonies(1).radiusMicron;
end

%restrict to colonies of the correct size
inds = [colonies.radiusMicron] == colSize;
colonies = colonies(inds);

chans = 1:length(meta.channelLabel);
chansToPlot = setdiff(chans,DAPIChannel);

r = imfilter(colonies(1).radialProfile.BinEdges,[1 1]/2)*meta.xres;
r(1) = colonies(1).radialProfile.BinEdges(1)*meta.xres;
r = r(1:end-1);
colCat = cat(3,colonies(:).radialProfile);

nucAvgAll = mean(cat(3,colCat.NucAvg),3);

nucAvgAllNormalized = bsxfun(@rdivide, nucAvgAll, nucAvgAll(:,DAPIChannel));

plot(r, nucAvgAllNormalized(:,chansToPlot),'.-','LineWidth',3)
legend(meta.channelLabel(chansToPlot));
axis([min(r) max(r) 0 max(max(nucAvgAllNormalized(:,chansToPlot)))]);

% how normalize cytoplasmic levels?%cytAvgAll = mean(cat(3,colCat.CytAvg),3);

%cytAvgAllNormalized = bsxfun(@rdivide, cytAvgAll, nucAvgAll(:,DAPIChannel));


% hold on
% %plot(r, cytAvgAllNormalized(:,2:4),'--')
% hold off