function plotInterColonyRadialVariability(meta, colonies, permutation, normalized)

radProfiles = cat(1,colonies.radialProfile);
radNucProfiles = cat(3, radProfiles.NucAvg);
r= radProfiles(1).BinEdges(1:end-1);

if ~exist('normalized','var') 
    normalized = false;
    n = 4;
    labels = meta.channelLabel;
end
if normalized
    n = 3;
    % make this general when DAPI is not the first channel
    radNucProfiles = radNucProfiles(:,2:4,:)./radNucProfiles(:,1,:);
    labels = meta.channelLabel(2:4);
end
if ~exist('permutation','var') || isempty(permutation)
    permutation = 1:n;
end

figure('Position',[0 0 1600 300]);
m = 1;

for i = 1:n
    subplot(m,n,i);
    plot(r, squeeze(radNucProfiles(:,permutation(i),:)))
    hold on
    plot(r, mean(squeeze(radNucProfiles(:,permutation(i),:)),2),'k','LineWidth',2)
    %plot(r, radialAvgLSM.nucAvg(:,i),'r','LineWidth',2)
    hold off
    title(labels(i));
    axis square
end

end