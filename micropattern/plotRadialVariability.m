function plotInterColonyRadialVariability(meta, colonies, permutation)

radProfiles = cat(1,colonies.radialProfile);
radNucProfiles = cat(3, radProfiles.NucAvg);
r= radProfiles(1).BinEdges(1:end-1);

figure('Position',[0 0 1600 300]);
m = 1;
n = 4;

if ~exist('permutation','var')
    permutation = 1:n;
end

for i = 1:n
    subplot(m,n,i);
    plot(r, squeeze(radNucProfiles(:,permutation(i),:)))
    hold on
    plot(r, mean(squeeze(radNucProfiles(:,permutation(i),:)),2),'k','LineWidth',2)
    %plot(r, radialAvgLSM.nucAvg(:,i),'r','LineWidth',2)
    hold off
    title(meta.channelLabel(i));
    axis square
end

end