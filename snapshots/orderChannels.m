function order = orderChannels(meta)

    order = [0 0 0 0];
    preferredOrder = {'CY5','RFP','GFP','DAPI'};
    
    for i = 1:4
        order(i) = find(strcmp(meta.channelNames,preferredOrder{i}));
    end
end