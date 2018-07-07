function [radialAvgNuc, r] = plotAveragesNoSegmentation(...
                        meta, colSize, DAPIChannel, colonies, doubleNormalize)

% doubleNormalize: boolean
% first normalize by DAPI, then scale all profiles from 0 to 1

if ~exist('doubleNormalize','var')
    doubleNormalize = false;
end

radialAvg = makeAveragesNoSegmentation(...
                meta, colSize, DAPIChannel, colonies);
                        
chans = 1:length(meta.channelLabel);
chansToPlot = setdiff(chans,DAPIChannel);

if doubleNormalize
    radialAvgNuc = radialAvg.nucAvgDAPImaxNormalized;
else
    radialAvgNuc = radialAvg.nucAvgDAPINormalized;
end
r = radialAvg.r;

plot(r, radialAvgNuc(:,chansToPlot),'.-','LineWidth',3)
axis([min(r) max(r) 0 max(max(radialAvgNuc(:,chansToPlot)))]);
legend(meta.channelLabel(chansToPlot));