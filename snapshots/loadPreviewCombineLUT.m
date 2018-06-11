function [preview, segpart, preview8bit, lim, Ilim] = loadPreviewCombineLUT(dataDir, filelist, dataChannels, nucChannel)

    n = 1;
    m = 1;
    I = ones([1 numel(filelist)]);
    ymin = I*n*2048;
    xmin = I*n*2048;
    ymax = ymin + m*2048;
    xmax = xmin + m*2048;
    preview = {};
    Ilim = {};
    segpart = {};

    series = 1;

    % load data and determine intensity limits per file
    for fi = 1:numel(filelist)

        vsifile = fullfile(dataDir,filelist{fi});
        metatmp = Metadata(vsifile);
        xmax(fi) = min(xmax(fi), metatmp.xSize); 
        ymax(fi) = min(ymax(fi), metatmp.ySize);
        W = xmax(fi)-xmin(fi)+1; H = ymax(fi)-ymin(fi)+1;
        img_bf = bfopen_mod(vsifile,xmin(fi),ymin(fi),W,H,series);

        time = 1;
        P = Position(numel(dataChannels), vsifile, time);
        P.setID(fi);
        seg = P.loadSegmentation(fullfile(dataDir,'MIP'), nucChannel);
        segpart{fi} = seg(ymin(fi):ymax(fi),xmin(fi):xmax(fi));
    
        for ci = 1:numel(dataChannels)

            im = img_bf{1}{dataChannels(ci),1};
            preview{ci, fi} = im;
            maxim = double(max(im(:)));
            minim = double(min(im(:)));
            tolerance = 0.04;
            Ilim{dataChannels(ci), fi} = stretchlim(mat2gray(im), tolerance)*(maxim-minim) + minim;
        end
    end
    
    % combine lookup tables and save RGB previews

    preview8bit = {}; 
    for fi = 1:numel(filelist)
        for ci = 1:numel(dataChannels)  

            Imin = min(min([Ilim{dataChannels(ci),:}]));
            Imax = max(max([Ilim{dataChannels(ci),:}]));
            preview8bit{ci,fi} = uint8((2^8-1)*mat2gray(preview{ci, fi}, [Imin Imax]));
        end
    end
    
    lim = struct('xmin',xmin,'xmax',xmax,'ymin',ymin,'ymax',ymax);

    % small combined preview is good enough:
    
%     for fi = 1:numel(filelist)
%         previewChannels = 1:3;
%         RGBim = cat(3,preview8bit{previewChannels,fi});
%         
%         [~,barefname,~] = fileparts(filelist{fi});
%         filename = fullfile(dataDir, 'preview', ['RGBpreview_' barefname '_' [meta.channelLabel{dataChannels}] '.tif']);
%         imwrite(RGBim, filename);
%     end    
end