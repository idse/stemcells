clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

% dataDir = '/Users/idse/data_tmp/160901_smad2';

mainDataDir = '/Users/idse/data_tmp/170207_hESCsmad2series/170207_siS4A/';
dataDirs = fullfile(mainDataDir, {...
    'A1','A2','A4','A6','A9','A12','A24','SB','X12','X24',...
    });
dataDir = dataDirs{1};

MIPdirs = fullfile(dataDirs,'MIP');

list = dir(fullfile(dataDirs{1},'*vsi'));
vsifile = fullfile(dataDirs{1},list(1).name);

meta = Metadata(vsifile);
%filenameFormat = meta.filename;

% determine process numbers of vsi files
processnrs = {};

for di = 1:numel(dataDirs)
    vsifiles = dir(fullfile(dataDirs{di},'*.vsi'));
    processNumbers = zeros([numel(vsifiles) 1],'uint16');
    for i = 1:numel(vsifiles)
        s = strsplit(vsifiles(i).name,{'_','.'});
        processNumbers(i) = uint16(str2double(s{2}));
    end
    processnrs{di} = processNumbers;
end

% manual metadata
%-------------------

meta.posPerCondition = 4;
meta.nWells = 10;
meta.nPositions = meta.nWells*meta.posPerCondition;

nucChannel = 1;
S2Channel = 2;
S4Channel = 3;

meta.channelLabel = {'H2B','Smad2','Smad4'};

%% extract prefixes for different conditions

% convention is conditionName_p0000.tif etc
% this extracts all the unique condition names in the dataDir

conditions = {};
for di = 1:numel(MIPdirs)
    
    listing = dir(fullfile(MIPdirs{di},'*tif'));
    
    for i = 1:numel(listing)

        filename = listing(i).name;
        [~, barefname] = fileparts(filename);
        s = strsplit(barefname,'_MIP');
        conditions = [conditions s{1}];
    end
end

meta.conditions = unique(conditions,'stable');

% %% previews
% 
% for di = 1:meta.nWells
%     
%     metatmp = meta;
%     metatmp.nWells = 1;
% 
%     stitchedPreviews(dataDirs{di}, metatmp);
% end

%% extract nuclear and cytoplasmic levels

% externally I will have all indices starting at 1
% the Andor offset to start at 0 will be internal

opts = struct(  'cytoplasmicLevels',    true,... %'tMax', 25,...
                    'dataChannels',     1:meta.nChannels,... 
                    'fgChannel',        S4Channel,...
                    'tMax',             meta.nTime,...
                    'nucShrinkage',     2,...
                    'cytoSize',         5,...
                    'cytoMargin',       5,...
                    'bgMargin',         4);

opts.cleanupOptions = struct('separateFused', true,'openSize',8,...
    'clearBorder',true, 'minAreaStd', 1, 'minSolidity',0, 'minArea',800);

%% finding the right options for extract data

di = 1;
dataDir = dataDirs{di};
opts.segmentationDir = fullfile(dataDir,'MIP');

pi = 1;
vsifile = fullfile(dataDir,['Process_' num2str(processnrs{di}(pi)) '.vsi']);
P = Position(meta.nChannels, vsifile, meta.nTime);
P.setID(pi);
seg = P.loadSegmentation(fullfile(dataDir,'MIP'), nucChannel);

% try out the nuclear cleanup settings on some frame:
time = 1;
bla = nuclearCleanup(seg(:,:,time), opts.cleanupOptions);
figure, imshow(cat(3,mat2gray(bla),seg(:,:,time),seg(:,:,time)))

%%

debugInfo = P.extractData(dataDir, nucChannel, opts);

%bgmask = debugInfo.bgmask;
nucmask = debugInfo.nucmask;
cytmask = false(size(nucmask));
cytmask(cat(1,debugInfo.cytCC.PixelIdxList{:}))=true;

bg = P.cellData(time).background;
nucl = P.cellData(time).nucLevelAvg;
cytl = P.cellData(time).cytLevelAvg;
(nucl-bg)/(cytl - bg)

im = P.loadImage(dataDir, S4Channel, time);
%im2 = P.loadImage(dataDir, nucChannel, time);

MIP = max(im,[],3);
A = imadjust(mat2gray(MIP));
s = 0.4;
imshow(cat(3, A + seg(:,:,time), A + s*nucmask, A + s*cytmask));

%%
tic
positions(meta.nPositions) = Position();
pi = 1;
for di = 1:meta.nWells
    for wpi = 1:meta.posPerCondition
    
        opts.segmentationDir = fullfile(dataDirs{di},'MIP');
        vsifile = fullfile(dataDirs{di},['Process_' num2str(processnrs{di}(wpi)) '.vsi']);

        positions(pi) = Position(meta.nChannels, vsifile, meta.nTime);
        positions(pi).setID(wpi);
        positions(pi).extractData(dataDirs{di}, nucChannel, opts);
        
        pi = pi + 1;
        
        save(fullfile(dataDir,'positions'), 'positions');
    end
end
toc

%%
load(fullfile(dataDir,'positions'));

%% plot smad2 vs time

figure,

ratioWellMeans = {};
ratioWell = {};

channels = 2:3;
normalize = false;
tocyt = false;

if tocyt
    ylabelstr = 'nuclear : cytoplasmic intensity';
else
    ylabelstr = 'nuclear intensity';
end    
if normalize
    ylabelstr = ['normalized ' ylabelstr];
end

for ci = channels
    
    pi = 1;
    ratioWellMeans{ci} = zeros([1 meta.nWells]);

    for wi = 1:meta.nWells

        ratioWell{ci,wi} = zeros([1 meta.posPerCondition]);

        for wpi = 1:meta.posPerCondition

            nucl = positions(pi).cellData.nucLevelAvg(ci);
            cytl = positions(pi).cellData.cytLevelAvg(ci);
            bg = positions(pi).cellData.background(ci);
            ratio = (nucl-bg);%/(cytl - bg);
            if tocyt
                ratio = ratio/(cytl - bg);
            end

            ratioWell{ci,wi}(wpi) = ratio;

            pi = pi + 1;
        end

        ratioWellMeans{ci}(wi) = mean(ratioWell{ci,wi});
    end
end

clf
hold on
colors = {'b','r'};
lw = 2;
labels = meta.channelLabel(opts.dataChannels);
times = [1 2 4 6 9 12 24 -1 -0.5 0];

for ci = 1:numel(channels)
    [~,order] = sort(times);
    ordertimes = times(order);
    
    rwm = ratioWellMeans{channels(ci)}(order);
    
    baseline = ratioWellMeans{channels(ci)}(times == -1);
    amplitude = max(rwm) - baseline;
    
    if normalize
        rwm = (rwm - baseline)./amplitude;
    end
    plot(ordertimes(2:end), rwm(2:end), '-', 'Color',colors{ci}, 'LineWidth',lw);
end

legend(labels(channels));

for ci = 1:numel(channels)
    for wi = 1:meta.nWells
        for wpi = 1:meta.posPerCondition
            rw = ratioWell{channels(ci),wi};
            
            baseline = ratioWellMeans{channels(ci)}(times == -1);
            amplitude = max(ratioWellMeans{channels(ci)}) - baseline;
            
            if normalize
                rw = (rw - baseline)./amplitude;
            end
            plot(times(wi), rw, '.', 'MarkerSize',10,...
                'Color',colors{ci}, 'LineWidth',lw);
        end
    end
end

fs = 15;
xlabel('time (hours)', 'FontSize',fs,'FontWeight','Bold');
ylabel(ylabelstr, 'FontSize',fs, 'FontWeight','Bold');
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
    
hold off

%export_fig(fullfile(dataDir,[ylabelstr '.pdf']),'-native -m2');
saveas(gcf, fullfile(dataDir,[ylabelstr '.png']));
saveas(gcf, fullfile(dataDir,[ylabelstr '.fig']));

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
