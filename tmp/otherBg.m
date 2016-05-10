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
for ti = 1%:tmax

    disp(['processing time ' num2str(ti)]);
    
    for pi = 1:meta.nPositions

        disp(['reading MIP ' num2str(pi)]);
        % gridSize 1 and 2 may be swapped, I have no way of knowing right now
        [i,j] = ind2sub(gridSize,pi);

        fname = fullfile(dataDir,'MIP',[barefname sprintf('_MIP_p%.4d_w%.4d.tif',pi-1,nucChannel-1)]);
        imgsNuc{j,i} = double(imread(fname,ti));

        fname = fullfile(dataDir,'MIP',[barefname sprintf('_MIP_p%.4d_w%.4d.tif',pi-1,S4Channel-1)]);
        imgsS4{j,i} = double(imread(fname,ti));
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

%%
Ilim = [500 5000];
A = imgsS4{3,2}(:,end-200+1:end);
B = imgsS4{3,3}(:,1:200);
[shifti,shiftj] = xcorr2fft(A,B);

Ap = A(1-shifti:end,1:end-shiftj);
Bp = B(1:end+shifti,shiftj+1:end);

imshow(cat(3,   mat2gray(Ap,Ilim),...
                mat2gray(Bp,Ilim),double(0*Ap)))
            
%%
size(A(:,1:end+shifti))
size(B(:,1-shifti:end))

%%
imshow(A,Ilim)
figure,
imshow(B,Ilim)

%%
plot(mean(Ap,2))
hold on
plot(mean(Bp,2),'-r')
hold off
%%
bla = medfilt2(mean(Ap,2)./mean(Bp,2),[5 1]);
plot(bla(1:end-1))

%%
plot(mean(Ap,1))
hold on
plot(mean(Bp,1),'-r')
plot(0*mean(Bp,1) + mean(mean(Bp,1)+mean(Ap,1),1)/2,'g')
hold off

%%
bg = 0;
bla = medfilt2(mean(Ap-bg,1)./mean(Bp-bg,1),[1 1]);
plot(bla(1:end-1))

%%
bg = 0;
App = mean(Ap,1) - (mean(Ap-bg,1)-mean(Bp-bg,1))/2;
Bpp = mean(Bp,1) + (mean(Ap-bg,1)-mean(Bp-bg,1))/2;
plot(mean(Ap,1), 'r')
hold on
plot(App, '--r')
plot(mean(Bp,1), 'g')
hold off

%%
App = bsxfun(@rdivide, double(Ap-bg), double(bla));
mean(App(:))/mean(Bp(:))
imshow(cat(3,mat2gray(App+bg,Ilim),mat2gray(Bp,Ilim),0*App),[])

%%

plot(bla)
x = 1:numel(bla);
c = (bla(end)-bla(1))/numel(bla);
hold on
plot(bla(end)-c*x(end) + c*x)
hold off

%%
plot(mean(Ap,1))
hold on
plot(mean(Bp,1),'-r')
scale = (bla(end)-c*x(1,end) + c*x/2);
scale = 1 + scale - mean(scale(:));
plot(mean(Ap,1)./scale,'c')
plot(mean(Bp,1).*scale,'m')
plot(0*mean(Bp,1) + mean(mean(Bp,1)+mean(Ap,1),1)/2,'g')
hold off

%%
%[X,Y] = meshgrid(1:1024,1:1024);
[X,Y] = meshgrid(1:size(Ap,2),1:size(Ap,1));
scale = (bla(end)-c*X(1,end) + c*X);
%scale = 1 + scale - mean(scale(:));
App = double(Ap-bg)./scale;
Bpp = double(Bp-bg);%./scale;
mean(App(:))/mean(Bpp(:))
imshow(cat(3,mat2gray(App+bg,Ilim),mat2gray(Bpp,Ilim),0*App),[])

%%
[X,Y] = meshgrid(1:1024,1:1024);
scale = (bla(end)-c*X(1,end) + c*X/2);
scale = 1 + scale - mean(scale(:));

Ax = imgsS4{1,1}./scale;
Ax = Ax(:,end-200+1:end);

Bx = imgsS4{1,2}./scale;
Bx = Bx(:,1:200);

Axp = Ax(1-shifti:end,1:end-shiftj);
Bxp = Bx(1:end+shifti,shiftj+1:end);

imshow(cat(3,   mat2gray(Axp,Ilim),...
                mat2gray(Bxp,Ilim),double(0*Ap)))

%%
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

%% extract nuclear and cytoplasmic levels

tic
for pi = 1%:meta.nPositions
    
    position = DynamicPositionAndor(meta, pi);
    
    opts = struct(  'cytoplasmicLevels',true,...
                    'tMax', 10,...
                    'dataChannels',     S4Channel,...
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'background',       ones([meta.xSize meta.ySize],'uint16'));

    position.extractData(dataDir, nucChannel, opts)
end
toc

%%

mean(position.cellData(1).nucLevel)

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

%save(fname, 's');

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

