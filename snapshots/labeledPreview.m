function labeledPreview(meta, preview8bit, previewChannels)

    N = numel(meta.conditions);
    n = ceil(N/4);
    %m = ceil(N/n);
    m = min(N,4);
    screensize = get( 0, 'Screensize' );
    margin = 50;
    fs = 15;
    w = screensize(3);
    h = n*(screensize(3)/m + margin/2);
    % if h > (screensize(4)-100)
    %     w = w*(screensize(4)-100)/h;
    %     h = screensize(4);
    % end
    figure('Position', [1, 1, w, h]);
    for fi = 1:numel(meta.conditions)

        subplot_tight(n,m,fi,[0.05 0])
        RGBim = cat(3,preview8bit{previewChannels,fi});
        %RGBnuc = repmat(preview8bit{dataChannels==nucChannel,fi},[1 1 3]);
        %RGBim = RGBim + 0.5*RGBnuc;
        imshow(RGBim);
        title(meta.conditions(fi),'FontSize',fs,'FontWeight','bold');
        labelstr = ['\color{red}'   meta.channelLabel{previewChannels(1)}...
                    '\color{green}' meta.channelLabel{previewChannels(2)}...
                    '\color[rgb]{0.1, 0.5, 1}' meta.channelLabel{previewChannels(3)}];
        text(margin, size(RGBim,1) - 2*margin, labelstr,'FontSize',fs,'FontWeight','bold');
    end
end

