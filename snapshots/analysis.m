clear all; close all;

[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(scriptPath)
addpath(scriptPath);

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

dataDir = '/Volumes/Seagate Backup Plus Drive/160317_ibidi_RIactivin';
meta = MetadataAndor(dataDir);
filenameFormat = meta.filename;
 
barefname = 'RIactivin100';
treatmentTime = 8; % first time point after treatment

nucChannel = 2;
S4Channel = 1;

% % visualize positions
% 
% meta.displayPositions;

%% read the MIPs from previous step at time 1 

gridSize = meta.montageGridSize;
pixelOverlap = round(1024*meta.montageOverlap/100);

imgsNuc = {};
imgsS4 = {};

tmax = meta.nTime;
for ti = 1:tmax

    disp(['processing time ' num2str(ti)]);
    
    for pi = 1:meta.nPositions

        disp(['reading MIP ' num2str(pi)]);
        % gridSize 1 and 2 may be swapped, I have no way of knowing right now
        [i,j] = ind2sub(gridSize,pi);

        fname = fullfile(dataDir,'MIP',[barefname sprintf('_MIP_p%.4d_w%.4d.tif',pi-1,nucChannel-1)]);
        imgsNuc{j,i} = imread(fname,ti);

        fname = fullfile(dataDir,'MIP',[barefname sprintf('_MIP_p%.4d_w%.4d.tif',pi-1,S4Channel-1)]);
        imgsS4{j,i} = imread(fname,ti);

    end
    
    % stitch together
    if ti == 1
        % get register positions of upper left corner
        upperleft = registerImageGrid(imgsNuc, pixelOverlap);
    end
    nucStitched = stitchImageGrid(upperleft, imgsNuc);
    S4Stitched = stitchImageGrid(upperleft, imgsS4);

    % make clean preview (not for quantitative analysis
    nucSmall = imfilter(nucStitched,[1 1]/2);
    nucSmall = nucSmall(1:2:end,1:2:end);
    nucSmall = imadjust(mat2gray(nucSmall));
    nucSmall = uint16((2^16-1)*nucSmall);

    S4Small = imfilter(S4Stitched,[1 1]/2);
    S4Small = S4Small(1:2:end,1:2:end);
    S4Small = imadjust(adapthisteq(mat2gray(medfilt2(S4Small,[3 3]))));
    S4Small = uint16((2^16-1)*S4Small);

    if ti == 1
        previewS4 = zeros([size(nucSmall) tmax],'uint16');
        previewNuc = zeros([size(nucSmall) tmax],'uint16');
    end
    previewNuc(:,:,ti) = nucSmall;
    previewS4(:,:,ti) = S4Small;
end

fname = fullfile(dataDir, 'stichedPreviewNuclei.tif');
if exist(fname,'file')
    delete(fname);
end
bfsave(previewNuc, fname);

fname = fullfile(dataDir, 'stichedPreviewS4.tif');
if exist(fname,'file')
    delete(fname);
end
bfsave(previewS4, fname);

% figure, imshow(cat(3,0*nucSmall,S4Small,0*S4Small));
% s = strsplit(meta.timeInterval,' ');
% dt = str2double(s{1});
% unit = s{2};
% t = (ti - treatmentTime)*dt;
% text(100,100,['T = ' num2str(t) s{2}],'Color','white','FontSize',18);

%% read nuclear segmentation

S4cyt = zeros([meta.nPositions meta.nTime]);
S4nuc = zeros([meta.nPositions meta.nTime]);

S4cytZmed = zeros([meta.nPositions meta.nTime]);
S4nucZmed = zeros([meta.nPositions meta.nTime]);
S4bgZmed = zeros([meta.nPositions meta.nTime]);

S4cytMIP = zeros([meta.nPositions meta.nTime]);
S4nucMIP = zeros([meta.nPositions meta.nTime]);
S4bgMIP = zeros([meta.nPositions meta.nTime]);

% pi and ti run from 1 here
for pi = 1%:meta.nPositions
    tic
    
    pi = 1;
    position = DynamicPositionAndor(meta, pi);
    
    seg = position.loadSegmentation(fullfile(dataDir,'MIP'), nucChannel);
    nucleiSeg = squeeze(seg(:,:,nucChannel,:));
    bgSeg = squeeze(seg(:,:,S4Channel,:));
    
    for ti = 1%:meta.nTime
        
        % read RFP MIP z-index
        nucleiMIPidx = position.loadMIPidx(fullfile(dataDir,'MIP'), nucChannel, ti);
        
        % read smad4 z-stack
        S4 = position.loadImage(dataDir, S4Channel, pi);
    
        % make masks
        %--------------------------------
        bgmask = ~imdilate(imclose(bgSeg(:,:,ti),strel('disk',10)),strel('disk',5));
        nucmaskraw = nucleiSeg(:,:,ti);
        options = struct('separateFused', true, 'clearBorder',true);
        nucmask = nuclearCleanup(nucmaskraw,options);

        % processing z-level per nucleus
        %--------------------------------
        CC = bwconncomp(nucmask);

        MIPidx = 0*nucleiMIPidx;
        for i = 1:numel(CC.PixelIdxList)
            MIPidx(CC.PixelIdxList{i}) = median(nucleiMIPidx(CC.PixelIdxList{i}));
        end
        %imshow(MIPidx,[])
 
        %%
        cytoIdx = imdilate(MIPidx, strel('disk',10));
        cytoIdx(nucmaskraw)=0;
        %imshow(cytoIdx,[])

        zimin = min(MIPidx(MIPidx>0));
        zimax = max(MIPidx(:));

        S4nucTmp = zeros([1 zimax]);
        S4cytTmp = zeros([1 zimax]);
        S4nucW = zeros([1 zimax]);
        S4cytW = zeros([1 zimax]);

        for zi = zimin:zimax

            submask = MIPidx == zi;
            if any(submask(:))
                S4slice = S4(:,:,zi);
                S4nucTmp(zi) = S4nucTmp(zi) + mean(S4slice(submask));
                S4nucW(zi) = S4nucW(zi) + sum(submask(:));

                submask = cytoIdx == zi;
                S4cytTmp(zi) = S4cytTmp(zi) + mean(S4slice(submask));
                S4cytW(zi) = S4cytW(zi) + sum(submask(:));
            end
        end

        S4nuc(pi,ti) = sum(S4nucTmp.*S4nucW)/sum(S4nucW);
        S4cyt(pi,ti) = sum(S4cytTmp.*S4cytW)/sum(S4cytW);
        
        % processing z-level median
        %--------------------------------
        zmed = median(MIPidx(MIPidx>0));

        S4slice = S4(:,:,zmed);
        %cytmask = imdilate(nucmask,strel('disk',10)) - nucmaskraw > 0;
        cytmask = cytoIdx > 0;
        S4nucZmed(pi,ti) = mean(S4slice(nucmask)); %S4nuc = mean(S4slice(MIPidx>0));
        S4cytZmed(pi,ti) = mean(S4slice(cytmask));
        S4bgZmed(pi,ti) = mean(S4slice(bgmask));
        
        % processing MIP
        %--------------------------------
        S4nucMIP(pi,ti) = mean(S4MIP(nucmask));
        S4cytMIP(pi,ti) = mean(S4MIP(cytmask));
        S4bgMIP(pi,ti) = mean(S4MIP(bgmask));
    end
    toc
end

%%

fname = fullfile(dataDir,'results.mat');

s = struct();
s.S4nucMIP = S4nucMIP;
s.S4cytMIP = S4cytMIP;
s.S4bgMIP = S4bgMIP;

s.S4nucZmed = S4nucZmed;
s.S4cytZmed = S4cytZmed;
s.S4bgZmed = S4bgZmed;

s.S4nuc = S4nuc;
s.S4cyt = S4cyt;

save(fname, 's');

%%

fname = fullfile(dataDir,'results.mat');
load(fname);

%%

t = (ones([1 16])'*(1:meta.nTime) - treatmentTime)'*10;
range = [min(t(1,:)) max(t(:,1)) 0.6 1.4];

figure(1)
N = s.S4nucMIP';
C = s.S4cytMIP';
bg = s.S4bgMIP';

plot(t,(N - bg)./(C - bg),'r')
axis(range)

figure(2)
N = s.S4nuc';
C = s.S4cyt';
bg = s.S4bgZmed';

plot(t,(N - bg)./(C - bg),'r')
axis(range)

figure(3)
N = s.S4nucZmed';
C = s.S4cytZmed';
bg = s.S4bgZmed';

plot(t,(N - bg)./(C - bg),'r')
axis(range)

%%

t = (ones([1 16])'*(1:meta.nTime) - treatmentTime)'*10;
range = [min(t(1,:)) max(t(:,1)) 0.6 1.4];

figure
N = s.S4nucMIP';
C = s.S4cytMIP';
bg = s.S4bgMIP';

plot(t,(N - bg)./(C - bg),'r')
axis(range)

hold on
N = s.S4nuc';
C = s.S4cyt';
bg = s.S4bgZmed';

plot(t,(N - bg)./(C - bg),'b')
axis(range)
hold off

%%

figure(4)
N = s.S4nucMIP';
C = s.S4cytMIP';
bg = s.S4bgMIP';

range = [min(t(1,:)) max(t(:,1)) min(s.S4bgZmed(:)) max(N(:))];

plot(t, N,'b')
hold on
plot(t,C,'g')
plot(t,bg,'c')
axis(range)
%plot(1000*(N - bg)./(C - bg),'r')
hold off

figure(5)
N = s.S4nuc';
C = s.S4cyt';
bg = s.S4bgZmed';

plot(t, N,'b')
hold on
plot(t,C,'g')
plot(t,bg,'c')
axis(range)
%plot(1000*(N - bg)./(C - bg),'r')
hold off

figure(6)
N = s.S4nucZmed';
C = s.S4cytZmed';
bg = s.S4bgZmed';

plot(t, N,'b')
hold on
plot(t,C,'g')
plot(t,bg,'c')
axis(range)
%plot(1000*(N - bg)./(C - bg),'r')
hold off

%% play with separating fused nuclei

% pi and ti run from 1 here
for pi = 1%:meta.nPositions
    
    fnameFormat = [barefname '_MIP_p%.4d_w%.4d_Simple Segmentation.h5'];
    fname = fullfile(dataDir, 'MIP', sprintf(fnameFormat,pi-1,nucChannel));
    nucleiSeg = squeeze(h5read(fname, '/exported_data')) == 2;
end

% make masks
%--------------------------------
ti = 4;
bgmask = ~imdilate(imclose(bgSeg(:,:,ti)',strel('disk',10)),strel('disk',5));
nucmaskraw = nucleiSeg(:,:,ti)';
nucmask = bwareaopen(nucmaskraw,500);
nucmask = imopen(nucmask,strel('disk',10));
        
imshow(nucmask)

% read RFP MIP
ci = nucChannel;
fnameFormat = [barefname '_MIP_p%.4d_w%.4d.tif'];
fname = fullfile(dataDir,'MIP',sprintf(fnameFormat,pi-1,ci));
nucleiMIP = imread(fname,ti);

%%
fti = ceil(ti/meta.tPerFile) - 1;
fname = fullfile(dataDir, sprintf(filenameFormat,pi-1,fti,S4Channel));
S4 = readStack(fname);
S4 = S4(:,:,:,:,meta.tPerFile-rem(ti,meta.tPerFile)); 
S4MIP = max(S4,[],3);

%%
newnucmask = separateFusedNuclei(nucmask);

%% segmentation overview

bgedge = bgmask - imerode(bgmask,strel('disk',3));
rgbim = label2rgb(bwlabel(newnucmask),'jet','k','shuffle');
S4comb = mat2gray(S4MIP);
rgbCombined = double(rgbim);
for i = 1:3
    nucOutline = rgbim(:,:,i) - imerode(rgbim(:,:,i),strel('disk',3));
    rgbCombined(:,:,i) = mat2gray(nucOutline) + bgedge + S4comb;
end
imshow(rgbCombined)

%% make distributions overall and binned by boundary distance
% start in single panel

