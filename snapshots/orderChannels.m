function order = orderChannels(meta)

    order = [0 0 0 0];
    preferredOrder = {'CY5','RFP','GFP','DAPI'};
    
    for i = 1:4
        j = find(strcmp(meta.channelNames,preferredOrder{i}));
        if ~isempty(j)
            order(i) = j;
        end
    end
end