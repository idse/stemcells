clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

%dataDir = '/Volumes/IdseData3/161009_Smad2FISH';
dataDir = '/Users/idse/data_tmp/161009_Smad2FISH';
MIPdir = fullfile(dataDir,'MIP');

meta = MetadataAndor(dataDir);
filenameFormat = meta.filename;

% manual metadata
%-------------------

meta.posPerCondition = 6;
meta.nWells = 7;
meta.nPositions = meta.nWells*meta.posPerCondition;

nucChannel = 4;
S4Channel = 1; % actually Smad2 here

meta.channelLabel = {'Smad2','Lefty (560)','Nodal (640)','DAPI'};

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
meta.conditions = meta.conditions([7 1 3 4 5 6 2]);

%% extract nuclear and cytoplasmic levels

% IS IT ACTUALLY LOADING THE WHOLE STACK OR DOES IT THINK THIS IS FLAT?

S2Channel = 4;
opts = struct(  'cytoplasmicLevels',    true,... %'tMax', 25,...
                    'dataChannels',     1:4,...
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'nucShrinkage',     2,...
                    'cytoSize',         5,...
                    'bgMargin',         4);
%                                     'fgChannel',        S4Channel,...
%                    'MIPidxDir',        fullfile(dataDir,'MIP'),...

%%
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
imshow(cat(3, A + 0*nucmask, A + s*nucmask, A + s*cytmask));

%%
load(fullfile(dataDir,'positions'));

%%
% FISH analysis
%----------------

seg = {};
FISHchannels = 2:3;
ratioWellMeans = {};
ratioWell = {};
spotmode = {};
spotmodemean = {};

for ci = 2%FISHchannels
    
    ti = 1;
    ratioWellMeans{ci} = zeros([1 meta.nWells]);
    spotmodemean{ci} = zeros([1 meta.nWells]);
    
    for wi = 1:meta.nWells

        ratioWell{ci,wi} = zeros([1 meta.posPerCondition]);
        spotmode{ci,wi} = zeros([1 meta.posPerCondition]);

        for wpi = 1:meta.posPerCondition
            
            pi = (wi-1)*meta.posPerCondition + wpi;
            seg = positions(pi).loadSegmentation(opts.segmentationDir, ci);
            
            im = positions(pi).loadImage(dataDir, ci, ti);
            im = max(im,[],3);
            im2 = im - imerode(im,strel('disk',10));
            
            ratioWell{ci,wi}(wpi) = sum(im2(seg))/positions(pi).ncells;
            
            %ratioWell{ci,wi}(wpi) = sum(seg(:))/positions(pi).ncells; 
            %ratioWell{ci,wi}(wpi) = (sum(im(seg)) - sum(seg(:))*mean(im(~seg)))/positions(pi).ncells;
            %ratioWell{ci,wi}(wpi) = (sum(im2(seg)) - sum(seg(:))*mean(im2(~seg)))/positions(pi).ncells;
 
            % from segmentation of spots, get all isolated spots by area
            % cutoff and determine the mode of the mean intensity distribution
            stats = regionprops(seg,'Area','PixelIdxList');
            stats = stats([stats.Area] < 9 & [stats.Area] > 4);
            bla = zeros([1 numel(stats)]);
            for i = 1:numel(stats)
                bla(i) = mean(im2(stats(i).PixelIdxList));
            end
            bins = 500:50:4000;
            n = hist(bla,bins);
            [~,idx] = max(n);
            spotmode{ci,wi}(wpi) = bins(idx);
        end

        ratioWellMeans{ci}(wi) = mean(ratioWell{ci,wi});
        spotmodemean{ci}(wi) = mean(spotmode{ci,wi});
    end
end

%% plot the spot modes and spot mode mean
% the expectation is that they are roughly the same at each time unless
% there are variation in imaging coniditons etc

figure, 
hold on
times = [-1 0 1 2 4 8 11];
colors = lines(4);

for ci = FISHchannels

    for wi = 1:meta.nWells
        plot(times(wi), spotmode{ci, wi},'.','Color',colors(ci,:))
        plot(times(wi), mean(spotmode{ci, wi}),'x','Color',colors(ci,:))
    end
end
%ylim([800 1300]);

%A*sqrt(pi/a)
%a = 1/2 s^2 -> A*sqrt(2 pi)*s
% res 320 nm / pixel? s ~ 1 pixel? or psf

%% plot smad2 vs time

for ci = 1
    
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

% %%
% for ci = 2:3%channels
%     
%     pi = 1;
%     ratioWellMeans{ci} = zeros([1 meta.nWells]);
% 
%     for wi = 1:meta.nWells
% 
%         ratioWell{ci,wi} = zeros([1 meta.posPerCondition]);
% 
%         for wpi = 1:meta.posPerCondition
% 
%             nucl = positions(pi).cellData.nucLevelAvg(ci);
%             cytl = positions(pi).cellData.cytLevelAvg(ci);
%             bg = positions(pi).cellData.background(ci);
%             ratio = (nucl + cytl - bg)/500;
% 
%             ratioWell{ci,wi}(wpi) = ratio;
% 
%             pi = pi + 1;
%         end
% 
%         ratioWellMeans{ci}(wi) = mean(ratioWell{ci,wi});
%     end
% end
%%
figure,
normalize = true;%false;

channels = 1:3;
clf
hold on
colors = {'b','r','g'};
labels = meta.channelLabel(channels);
times = [-1 0 1 2 4 8 11];

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

%% just the FISH data 
figure,
normalize = true;%false;

channels = 2:3;

colors = {'b','r','g'};
labels = meta.channelLabel(channels);
times = [-1 0 1 2 4 8 11];
N = 7;

for ci = channels
    [~,order] = sort(times);
    rwm = ratioWellMeans{ci}(order);
    if normalize
        %rwm = rwm./rwm(2);% normalize by time 0
        rwm = rwm./(spotmodemean{ci}*N); % normalize by single spot intensity
    end
    plot(times(order), rwm, '-', 'Color',colors{ci});
    hold on
end

legend(labels, 'Location','SouthEast');

for ci = channels
    for wi = 1:meta.nWells
        for wpi = 1:meta.posPerCondition
            rw = ratioWell{ci,wi};
            if normalize
                %rw = rw./ratioWellMeans{ci}(2); % normalize by time 0
                rw = rw./(spotmodemean{ci}(wi)*N); % normalize by single spot intensity
            end
            plot(times(wi), rw, '.', 'Color',colors{ci});
        end
    end
end

if normalize
    ylabelstr = 'approximate spot count per cell';
else
    ylabelstr = 'probe intensity per cell';
end

fs = 15;
xlabel('time (hours)', 'FontSize',fs,'FontWeight','Bold');
ylabel(ylabelstr, 'FontSize',fs, 'FontWeight','Bold');
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
    
hold off


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
