clear all; close all;

clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
dataDir = '/Users/idse/data_tmp/0_kinetics/170519_FRAPctrl';

oibfiles = {'170518_FrapContr2.oif','170518_FrapContr3.oif',...
            '170518_FrapContr4.oif'};

tmaxall = {[140 60],[],[32 32]};

allresults = {};

descriptions = {'untreated', 'untreated', 'untreated'};

nucChannel = 2;
S4Channel = 1;

%%

for fi = 1:numel(oibfiles)
    
    oibfile = oibfiles{fi};
    [~,barefname,~] = fileparts(oibfile);
    barefname = strrep(barefname,'.','dot');
    
    r = bfGetReader(fullfile(dataDir,oibfile));
    omeMeta = r.getMetadataStore();

    tcut = r.getSizeT()-1; % cutoff time
    tres = double(omeMeta.getPixelsTimeIncrement(0).value);
    Nfrapped = omeMeta.getShapeCount(0)/2;
    if isempty(tmaxall{fi})
        tmaxall{fi} = tcut*ones([1 Nfrapped]);
    end

    % read the data
    channels = 1:2;
    img = zeros([r.getSizeY() r.getSizeX() numel(channels) r.getSizeZ() tcut], 'uint16');
    for ti = 1:tcut
        for cii = 1:numel(channels)
            for zi = 1:r.getSizeZ()
                img(:,:,cii,zi,ti) = bfGetPlane(r, r.getIndex(zi-1,channels(cii)-1,ti-1)+1);
            end
        end
    end
    r.close();

    % set up results struct
    results = struct('description', descriptions{fi}, 'tmax', tmaxall(fi));
    results.tres = tres;
    results.xyres = double(omeMeta.getPixelsPhysicalSizeX(0).value);
    results.barefname = barefname;
    
    % read FRAP regions
    [x,y,shapeTypes] = readFRAPregions(omeMeta);
    results.shapeTypes = shapeTypes;

    % read out profile in the mask
    data = squeeze(img(:,:,S4Channel,zi,:));
    [tracesNuc, tracesNucNorm, maskall] = readFRAPprofile(data,x,y);
    results.tracesNuc = tracesNuc;
    results.tracesNucNorm = tracesNucNorm;
    
    % fitting
    results.bleachType = 'nuclear';
    results.fitType = 'nuclear';
    results.tmax = tmaxall{fi};
    
    [parameters, frapframe, gof] = fitFRAP(results);
    results.A = parameters.A;
    results.k = parameters.k;
    results.frapframe = frapframe;
    
    % store results of this video
    allresults{fi} = results;
    
    % visualize FRAP regions for diagnostic purposes

    % initiale state
    colors = lines(Nfrapped);
    ti = 1;
    zi = 1;
    ci = S4Channel; 
    imshow(imadjust(mat2gray(img(:,:,ci,zi,ti))))
    hold on
    for shapeIdx = 1:Nfrapped
        plot(x{shapeIdx}, y{shapeIdx},'Color',colors(shapeIdx,:),'LineWidth',2)
    end
    hold off
    saveas(gcf, fullfile(dataDir, [barefname '_FRAPframe' num2str(ti) '.png']));
    close;

    % final state
    ti = tcut;
    zi = 1; 
    IS4 = imadjust(mat2gray(img(:,:,ci,zi,ti)));
    Inuc = imadjust(mat2gray(img(:,:,nucChannel,zi,ti)));
    figure, imshow(cat(3, IS4 + Inuc, IS4 + maskall - imerode(maskall,strel('disk',3)), IS4));
    saveas(gcf, fullfile(dataDir, [barefname '_FRAPframeFinal.png']));
    close;
end

save(fullfile(dataDir,'results'),'allresults');

%% visualize results

for fi = 1:numel(oibfiles)
    
    tracesnorm = allresults{fi}.tracesNucNorm;
    tcut = size(tracesnorm,2);
    tres = allresults{fi}.tres;
    
    % filename
    oibfile = oibfiles{fi};
    [~,barefname,~] = fileparts(oibfile);
    barefname = strrep(barefname,'.','dot');
    
    % FRAP curves
    figure,
    t = repmat((1:tcut)*tres,[Nfrapped 1]);
    plot(t' ,tracesnorm');
    xlabel('time (sec)');
    ylabel('intensity')
    saveas(gcf,fullfile(dataDir, ['FRAPcurvesRaw_' barefname]));
    saveas(gcf,fullfile(dataDir, ['FRAPcurvesRaw_' barefname '.png']));

    plot(t' ,tracesnorm');
    xlabel('time (sec)');
    ylabel('normalized intensity')
    saveas(gcf,fullfile(dataDir, ['FRAPcurvesNorm_' barefname]));
    saveas(gcf,fullfile(dataDir, ['FRAPcurvesNorm_' barefname '.png']));
    close;

    % FRAP fit
    visualizeFRAPfit(allresults{fi})
    saveas(gcf,fullfile(dataDir, ['FRAPfit_' barefname]));
    saveas(gcf,fullfile(dataDir, ['FRAPfit_' barefname '.png']));
    close;
end


