%--------------------------------
% GENERATE OVERVIEW OF 24-WELL
%--------------------------------

clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

dataDir = '/Users/Idse/data_tmp/161217_siRNAoptimization';

fnameFormat = 'Process_%d.vsi';
fi0 = 1252;

vsifile = fullfile(dataDir, sprintf(fnameFormat, fi0 + 1));
[~,barefname,~] = fileparts(vsifile);

%%
% metadata
%----------------

meta = Metadata(vsifile);

meta.channelLabel = {'H2B','Smad4'};
nucChannel = 1;

%% assemble images from 24 well

% lay-out: 
% column wise down, 2 pics/well, 3 cols
picsPerWell = 2; 
wellsPerCol = 4;
picsPerCol = wellsPerCol*picsPerWell;

n = 1;
m = 1;
ymin = 1;
xmin = 1;
ymax = ymin + m*2048 - 1;
xmax = xmin + m*2048 - 1;

imgs = {};

tic

for fi = 1:24
    
    series = 1;
    ci = 1;
    vsifile = fullfile(dataDir, sprintf(fnameFormat, fi0 + fi));
    
    col = floor((fi-1)/picsPerCol);
    imi = mod(fi,2);
    row = rem(ceil(fi/picsPerWell)-1,wellsPerCol)+1;
    
    %[row, 2*col + 1 + imi]
    
    img_bf = bfopen_mod(vsifile,xmin,ymin,xmax-xmin+1,ymax-ymin+1,series,ci);
    im = img_bf{1}{1,1};
    
    % binning
    idx1 = 1:2:2048;
    idx2 = 2:2:2048;
    imbin = im(idx1,idx1) + im(idx1,idx2) + im(idx2,idx1) + im(idx2,idx2);
    
    imgs{row, 2*col + 1 + imi} = imbin;
end

toc

%%
imcat = {};
for i = 1:size(imgs,1)
    imcat{i} = cat(2,imgs{i,1:end});
end
imcat = cat(1,imcat{:});
%%
figure,
imcatAdjusted = mat2gray(imcat,[2000 8000]);
imshow(imcatAdjusted)

%%
imwrite(imcatAdjusted, fullfile(dataDir,'combined.tif'));

%%
imgs = {};

for fi = 25:26
    
    
    series = 1;
    ci = 1;
    vsifile = fullfile(dataDir, sprintf(fnameFormat, fi0 + fi))
    
    col = floor((fi-1)/picsPerCol);
    imi = mod(fi,2);
    row = rem(ceil(fi/picsPerWell)-1,wellsPerCol)+1;
    
    %[row, 2*col + 1 + imi]
    
    img_bf = bfopen_mod(vsifile,xmin,ymin,xmax-xmin+1,ymax-ymin+1,series,ci);
    im = img_bf{1}{1,1};
    
    % binning
    idx1 = 1:2:2048;
    idx2 = 2:2:2048;
    imbin = im(idx1,idx1) + im(idx1,idx2) + im(idx2,idx1) + im(idx2,idx2);
    
    imgs{1, 1 + imi} = imbin;
    
end

ctrlim = cat(2,imgs{1,1:end});
ctrlim = mat2gray(ctrlim,[2000 8000]);

imwrite(ctrlim, fullfile(dataDir,'ctrl.tif'));