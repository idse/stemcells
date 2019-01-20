clear all; close all;

addpath(genpath('/Users/idseimac/stemcells')); 
addpath(genpath('/Users/idse/repos/Warmflash/stemcells'));

dataDir = '/data/BCAT/D4';

meta = MetadataMicropatternAndor(dataDir);

filenameFormat = meta.filename;

% manual metadata
%-------------------

% TODO : modify MetadataAndor to contain all info below

treatmentTime = 2;
% workaround to not have to deal with changing stitched previews, actually
% it is 4x4
meta.posPerCondition = 4; 
meta.nWells = 16;
dataidx = {2:4, 5:7, 9:12, 13:16};
conditions = {'SB','SB_BMP50_iwp2','SB_BMP50','SB_BMP1'};

% dataidx = dataidx([1 4]);
% conditions = conditions([1 4]);

nucChannel = 2;
S4Channel = 1;

meta.montageGridSize = [2 2];
meta.montageOverlap = 57; % percentage 

tmax = meta.nTime-1;

% visualize positions
%---------------------

% meta.displayPositions;

% TODO: create merged cellData for montage
% movies of distribution over time

%% stitching

% meta2 = meta;
% meta2.nTime = meta.nTime - 1;
% ss = 1; % fullsize previews
% margin = 100;
% colonies = 1:24;%:12;%1:meta.nWells;
% type = 'MIP';
% saveStitchFullSize = true;
% upperleft = stitchedPreviews(dataDir, meta2, type, colonies, margin, ss, saveStitchFullSize); 
% 
% %save(fullfile(dataDir,'upperleft'),'upperleft');
% 
% %upperleft = stitchedPreviews(dataDir, meta); 
% 
% %%
% 
% meta2 = meta;
% meta2.nTime = meta.nTime - 1;
% colnr = 1;%1:meta.nWells;
% separatez = true;
% margin = 100;
% load(fullfile(dataDir,'upperleft'),'upperleft');
% 
% for colnr = [21:24 9:20]
%     makeStitchedStack(dataDir, meta2, upperleft, colnr, margin, separatez);
% end

%% create radial profiles using dirty segmentation

colfile = fullfile(dataDir, 'colonies.mat');
if exist(colfile,'file')
    disp('loading colonies');
    load(colfile, 'colonies');
else
    colonies = Colony;
end

% metadata
meta.colRadiiMicron = 350;
meta.colRadiiPixel = round(meta.colRadiiMicron/meta.xres);
dataParam = struct( 'dataType',         'BCAT',...
                    'meta',             meta,...
                    'dataDir',          dataDir,...
                    'fnameformat',      'stitched_well%d_w0001.tif',...
                    'segfnameformat',    'stitched_well%d_w0001_Simple segmentation.h5',...
                    'tmax',             tmax,...
                    'colonies',         [dataidx{:}]);
                
% processing parameters
findColParam = struct('sclose', 6, 'sopen', 8, 'checkcontained', false,...
                            'minArea', [],'convhull', true);
segMaskParam = struct('memDilate', 1);

colonies = dirtySegLive(dataParam, findColParam, segMaskParam, colonies);
save(fullfile(dataDir, 'colonies.mat'), 'colonies');

%% kymograph: 

kymographs = {};
param = struct('label','beta-catenin signaling','Irange',[0.6 0.8],...
                'xres',meta.xres,'tres',30, 'revbins', false);
            
for coli = [dataidx{:}]

    if ~isempty(colonies(coli).radiusPixel)
        
        makeKymograph(colonies(coli), param);
        saveas(gcf, sprintf(fullfile(dataDir,'kymographNorm_col%d.png'),coli));
        close;
    end
end

% make averages
for condi = 1:numel(conditions)

    makeKymograph(colonies(dataidx{condi}), param);
    saveas(gcf, fullfile(dataDir,['kymographNorm_' conditions{condi} '.png']));
    close;
end

%% plot radial profiles at different times as graphs

condi = 1;
[kymograph, ybins] = makeKymograph(colonies(dataidx{condi}), param);

clf
hold on
colors = jet(tmax);
for ti = 1:5:tmax
    plot(ybins, kymograph(ti,:),'Color',colors(ti,:),'LineWidth',2);
end
xlim([min(bins) max(bins)])
hold off

