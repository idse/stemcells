function batchMIP_Andor(inputdir, outputdir, channels, saveidx)

    meta = MetadataAndor(inputdir);
    
    for ci = 1:numel(channels)
        for pi = 1:meta.nPositions
            disp(['processing position ' num2str(pi) ', channel ' num2str(channels(ci))]);
            makeMIP_Andor(inputdir, pi-1, channels(ci), outputdir, saveidx(ci));
        end
    end
end