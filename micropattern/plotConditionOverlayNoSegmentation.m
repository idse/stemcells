function [radialAvgNuc, r] = plotConditionOverlayNoSegmentation(...
                        meta, colSize, DAPIChannel, colonies, conditions, doubleNormalize)

% doubleNormalize: boolean
% first normalize by DAPI, then scale all profiles from 0 to 1 on the same
% scale

% PAY ATTENTION
% colonies here is a cell array of arrays, with the cells corresponding to
% conditions

if ~exist('doubleNormalize','var')
    doubleNormalize = false;
end

nConditions = numel(colonies);
radialAvgNuc = {};
r = {};
minI = Inf*(1:meta.nChannels);
maxI = 0*(1:meta.nChannels);

for i = 1:nConditions
    
    radialAvg = makeAveragesNoSegmentation(...
                    meta, colSize, DAPIChannel, colonies{i});

    chans = 1:length(meta.channelLabel);
    chansToPlot = setdiff(chans,DAPIChannel);
    
    if ~isempty(DAPIChannel)
        radialAvgNuc{i} = radialAvg.nucAvgDAPINormalized;
    else
        radialAvgNuc{i} = radialAvg.nucAvg;
    end
    r{i} = radialAvg.r;    
    
    % for overall normalization
    % throw out 2 bins from edge when setting LUT
    % to prevent setting minimum by areas without cells
    Imargin = 6; 
    minI = min(minI, min(radialAvgNuc{i}(1:end-Imargin,:)));
    maxI = max(maxI, max(radialAvgNuc{i}(1:end-Imargin,:)));
end

if doubleNormalize
    for i = 1:nConditions
        for ci = 1:meta.nChannels
            radialAvgNuc{i}(:,ci) = (radialAvgNuc{i}(:,ci) - minI(ci))/(maxI(ci)-minI(ci));
        end
    end
end

colors = winter(nConditions);

m = 1;
for i = 1:meta.nChannels
    
    subplot_tight(m,meta.nChannels,i,0.02)
    hold on
    for j = 1:nConditions
        plot(r{j}, radialAvgNuc{j}(:,i),'.-','LineWidth',3,'Color',colors(j,:))
    end
    hold off
    axis([min(r{j}) max(r{j}) 0 1]);
    legend(conditions,'location','southwest');
    title(meta.channelLabel(i))
    
    axis square
%     if i > 1
%         legend off;
%     end
end

end