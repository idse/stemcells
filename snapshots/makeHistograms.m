function [bins, allHistograms] = makeHistograms(nucLevel, lim)

    allHistograms = {};
    bins = {};
    for fi = 1:numel(nucLevel)
        for channelIndex = 1:numel(lim)

            bins{channelIndex} = linspace(lim{channelIndex}(1),lim{channelIndex}(2),40);
            n = histc(nucLevel{fi}(:,channelIndex), bins{channelIndex});
            allHistograms{fi, channelIndex} = n;
        end
    end
end