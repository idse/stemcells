clear all; close all;

clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
dataDir = '/Users/idse/data_tmp/170209_FRAP';

oibfile = 'FRAP_A100nuclei1.5h.oib';

nucChannel = 2;
S4Channel = 1;

r = bfGetReader(fullfile(dataDir,oibfile));
omeMeta = r.getMetadataStore();

%% read the data

channels = 1:2;

img = zeros([r.getSizeY() r.getSizeX() numel(channels) r.getSizeZ() r.getSizeT()], 'uint16');
for ti = 1:r.getSizeT();
    for cii = 1:numel(channels)
        for zi = 1:r.getSizeZ()
            img(:,:,cii,zi,ti) = bfGetPlane(r, r.getIndex(zi-1,channels(cii)-1,ti-1)+1);
        end
    end
end

%% read a FRAP region and visualize

ROIidx = 0; % omeMeta.getROICount() -> only one ROI for multiple FRAP
shapeIdx = 0; % but multiple shapes

s = strsplit(char(omeMeta.getPolygonPoints(ROIidx,shapeIdx)),{' ',','});
s = cellfun(@str2double, s,'UniformOutput',false);
s = cat(1, s{:});

x = s(1:2:end);
y = s(2:2:end);

%% visualize

ti = 100;
zi = 1;
ci = S4Channel; %nucChannel;
imshow(img(:,:,ci,zi,ti),[1 2^12])
hold on
plot(x,y,'-g')
hold off

%% read out profile in the mask

mask = poly2mask(x, y, r.getSizeX(), r.getSizeY());

ci = S4Channel;

val = zeros([1 r.getSizeT()]);
for ti = 1:r.getSizeT()
    im = img(:,:,ci,zi,ti);
    val(ti) = mean(im(mask));
end

plot(val)

%%

r.close();
