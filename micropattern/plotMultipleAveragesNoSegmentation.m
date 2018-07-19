function [radialAvgNuc, r] = plotMultipleAveragesNoSegmentation(...
                        meta, colSize, DAPIChannel, colonies, conditions, doubleNormalize)

% doubleNormalize: boolean
% first normalize by DAPI, then scale all profiles from 0 to 1 on the same
% scale

if ~exist('doubleNormalize','var')
    doubleNormalize = false;
end

m = 1;
n = numel(colonies);
radialAvgNuc = {};
r = {};
minI = Inf*(1:meta.nChannels);
maxI = 0*(1:meta.nChannels);

for i = 1:n
    
    radialAvg = makeAveragesNoSegmentation(...
                    meta, colSize, DAPIChannel, colonies{i});

    chans = 1:length(meta.channelLabel);
    chansToPlot = setdiff(chans,DAPIChannel);

    radialAvgNuc{i} = radialAvg.nucAvgDAPINormalized;
    r{i} = radialAvg.r;    
    
    % for overall normalization
    minI = min(minI, min(radialAvg.nucAvgDAPINormalized));
    maxI = max(maxI, max(radialAvg.nucAvgDAPINormalized));
end

if doubleNormalize
    for i = 1:n
        for ci = 1:meta.nChannels
            radialAvgNuc{i}(:,ci) = (radialAvgNuc{i}(:,ci) - minI(ci))/(maxI(ci)-minI(ci));
        end
    end
end
            
for i = 1:n
    
    subplot_tight(m,n,i,0.02)
    plot(r{i}, radialAvgNuc{i}(:,chansToPlot),'.-','LineWidth',3)
    axis([min(r{i}) max(r{i}) 0 1]);
    legend(meta.channelLabel(chansToPlot));
    title(conditions{i})
    
    axis square
    if i > 1
        legend off;
    end
end

end