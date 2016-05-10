function chunkStruct = makeChunks(meta, maxMemoryGB, overlapWidth)
    % split the image up in big chunks 
    %
    % chunkStruct = makeChunks(meta, maxMemoryGB, overlapWidth)
    %
    % TODO: not hardcode bytesPerPixel fo 16 bit

    chunkStruct = struct();

    maxBytes = (maxMemoryGB*1024^3);
    
    bytesPerPixel = 2; 
    dataSize = meta.ySize*meta.xSize*meta.nChannels*bytesPerPixel;
    nChunks = ceil(dataSize/maxBytes);

    if nChunks > 1
        nRows = 2;
    else
        nRows = 1;
    end
    nCols = ceil(nChunks/nRows);

    xedge = (0:nCols)*(meta.xSize/nCols);
    yedge = (0:nRows)*(meta.ySize/nRows);

    xlim = {};
    ylim = {};
    
    for n = 1:numel(yedge)-1
        for m = 1:numel(xedge)-1
            
            xlim{n,m} = [uint32(xedge(m) + 1), uint32(xedge(m+1))];
            ylim{n,m} = [uint32(yedge(n) + 1), uint32(yedge(n+1))];
            
            if n < nRows 
                ylim{n,m}(2) = ylim{n,m}(2) + overlapWidth; 
            end
            if m < nCols
                xlim{n,m}(2) = xlim{n,m}(2) + overlapWidth;
            end
            
            height{n,m} = ylim{n,m}(2) - ylim{n,m}(1) + 1;
            width{n,m} = xlim{n,m}(2) - xlim{n,m}(1) + 1;
        end
    end
    
    chunkStruct.nChunks = nChunks;
    chunkStruct.nRows = nRows;
    chunkStruct.nCols = nCols;
%     chunkStruct.xedge = xedge;
%     chunkStruct.yedge = yedge;
    chunkStruct.xlim = xlim;
    chunkStruct.ylim = ylim;
    chunkStruct.width = width;
    chunkStruct.height = height;
end