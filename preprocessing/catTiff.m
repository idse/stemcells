function catTiff(infiles,outfile)

for fi = 1:numel(infiles)
    reader = bfGetReader(infiles{fi});
    %nT = reader.getSizeT;
    for ti = 1:numel(imfinfo(infiles{fi}))%nT
        %iplane = reader.getIndex(0,0,ti - 1) + 1;
        iplane = ti;
        img = bfGetPlane(reader,iplane);
        if ti==1 && fi==1
            imwrite(img,outfile,'Compression','none');
        else
            imwrite(img,outfile,'writemode','append','Compression','none');
        end
    end
end

end