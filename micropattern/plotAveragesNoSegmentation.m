function [radialAvgNuc, r] = plotAveragesNoSegmentation(...
                        datadir, colSize, DAPIChannel, filenr, doubleNormalize)

% doubleNormalize: boolean
% first normalize by DAPI, then scale all profiles from 0 to 1

[radialAvgNuc, r] = makeAveragesNoSegmentation(...
                datadir, colSize, DAPIChannel, filenr, doubleNormalize);
            
load(fullfile(datadir,'metaData.mat'));
            
chans = 1:length(meta.channelLabel);
chansToPlot = setdiff(chans,DAPIChannel);

plot(r, radialAvgNuc(:,chansToPlot),'.-','LineWidth',3)
axis([min(r) max(r) 0 max(max(radialAvgNuc(:,chansToPlot)))]);
legend(meta.channelLabel(chansToPlot));