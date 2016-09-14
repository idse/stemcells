clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

dataDir = '/Users/idse/data_tmp/160901_smad2';
MIPdir = fullfile(dataDir,'MIP');

meta = MetadataAndor(dataDir);
filenameFormat = meta.filename;

% manual metadata
%-------------------

meta.posPerCondition = 7;
meta.nWells = 7;
meta.nPositions = meta.nWells*meta.posPerCondition;

nucChannel = 1;
S4Channel = 2;

meta.channelLabel = {'DAPI','Smad4','H2B', 'Smad2'};

%% extract prefixes for different conditions

% convention is conditionName_p0000.tif etc
% this extracts all the unique condition names in the dataDir

listing = dir(fullfile(dataDir,'*tif'));

conditions = {};
for i = 1:numel(listing)
    
    filename = listing(i).name;
    [~, barefname] = fileparts(filename);
    s = strsplit(barefname,'_p');
    conditions = [conditions s{1}];
end

meta.conditions = sort(unique(conditions));

%% extract nuclear and cytoplasmic levels

% IS IT ACTUALLY LOADING THE WHOLE STACK OR DOES IT THINK THIS IS FLAT?

S2Channel = 4;
opts = struct(  'cytoplasmicLevels',    true,... %'tMax', 25,...
                    'dataChannels',     [S4Channel S2Channel],...
                    'fgChannel',        S4Channel,...
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'MIPidxDir',        fullfile(dataDir,'MIP'),...
                    'nucShrinkage',     2,...
                    'cytoSize',         5,...
                    'bgMargin',         4);

positions(meta.nPositions) = DynamicPositionAndor();
pi = 1;
for wi = 1:meta.nWells
    for wpi = 1:meta.posPerCondition
        meta.filename = [meta.conditions{wi} '_p%.4d.tif'];
        positions(pi) = DynamicPositionAndor(meta, wpi+1);
        positions(pi).extractData(dataDir, nucChannel, opts);
        pi = pi + 1;
    end
end

save(fullfile(dataDir,'positions'), 'positions');

%% finding the right options for extract data

wi = 1;
pi = 1;
meta.filename = [conditions{wi} '_p%.4d.tif'];
P = DynamicPositionAndor(meta, pi+1);
time = 1;
opts.tMax = time;

% try out the nuclear cleanup settings on some frame:
% bla = nuclearCleanup(seg(:,:,time), opts.cleanupOptions);
% imshow(bla)

debugInfo = P.extractData(dataDir, nucChannel, opts);

bgmask = debugInfo.bgmask;
nucmask = debugInfo.nucmask;
cytmask = false(size(nucmask));
cytmask(cat(1,debugInfo.cytCC.PixelIdxList{:}))=true;

bg = P.cellData(time).background
nucl = P.cellData(time).nucLevelAvg
cytl = P.cellData(time).cytLevelAvg
(nucl-bg)/(cytl - bg)

im = P.loadImage(dataDir, S4Channel, time);
MIP = max(im,[],3);
A = imadjust(mat2gray(MIP));
s = 0.4;
imshow(cat(3, A + 0*bgmask, A + s*nucmask, A + s*cytmask));

%%
load(fullfile(dataDir,'positions'));

%% plot smad2 vs time

ratioWellMeans = {};
ratioWell = {};

channels = 1:2;
normalize = false;

for ci = channels
    
    pi = 1;
    ratioWellMeans{ci} = zeros([1 meta.nWells]);

    for wi = 1:meta.nWells

        ratioWell{ci,wi} = zeros([1 meta.posPerCondition]);

        for wpi = 1:meta.posPerCondition

            nucl = positions(pi).cellData.nucLevelAvg(ci);
            cytl = positions(pi).cellData.cytLevelAvg(ci);
            bg = positions(pi).cellData.background(ci);
            ratio = (nucl-bg)/(cytl - bg);

            ratioWell{ci,wi}(wpi) = ratio;

            pi = pi + 1;
        end

        ratioWellMeans{ci}(wi) = mean(ratioWell{ci,wi});
    end
end

clf
hold on
colors = {'b','r'};
labels = meta.channelLabel(opts.dataChannels);
times = [0 -1 1 10.5 2 4 8];

for ci = channels
    [~,order] = sort(times);
    rwm = ratioWellMeans{ci}(order);
    if normalize
        rwm = rwm./max(rwm);
    end
    plot(times(order), rwm, '-', 'Color',colors{ci});
end

legend(labels(channels));

for ci = channels
    for wi = 1:meta.nWells
        for wpi = 1:meta.posPerCondition
            rw = ratioWell{ci,wi};
            if normalize
                rw = rw./max(ratioWellMeans{ci});
            end
            plot(times(wi), rw, '.', 'Color',colors{ci});
        end
    end
end

fs = 15;
xlabel('time (hours)', 'FontSize',fs,'FontWeight','Bold');
ylabel('nuclear : cytoplasmic intensity', 'FontSize',fs, 'FontWeight','Bold');
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
    
hold off

%export_fig(fullfile(dataDir,'Smad4vsSmad2.pdf'),'-native -m2');


%% distributions at each time

channels = 1:2;

for ci = channels
    
    pi = 1;

    for wi = 1:meta.nWells

        ratioWell{ci,wi} = [];

        for wpi = 1:meta.posPerCondition

            nucl = positions(pi).cellData.nucLevel(:,ci);
            cytl = positions(pi).cellData.cytLevel(:,ci);
            bg = positions(pi).cellData.background(ci);
            ratio = (nucl-bg)./(cytl - bg);

            ratioWell{ci,wi} = cat(1,ratioWell{ci,wi}, ratio);

            pi = pi + 1;
        end
    end
end

%%
ci = 2;
%x = 0.2:0.1:2.5;   % Smad4
x = 0.2:0.2:8;      % Smad2
colors = jet(meta.nWells);
lw = 2;
clf 
hold on
for wi = 1:meta.nWells
    n = hist(ratioWell{ci,order(wi)}, x);
    n = n./sum(n);
    plot(x,n,'Color',colors(wi,:),'LineWidth',lw)   
end
hold on
legend(meta.conditions(order))
xlim([x(1) x(end)]);

fs = 15;
xlabel(['nuclear : cytoplasmic ' labels{ci}], 'FontSize',fs,'FontWeight','Bold');
ylabel('frequency', 'FontSize',fs, 'FontWeight','Bold');
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')

export_fig(fullfile(dataDir,'Smad2distributions.pdf'),'-native -m2');
