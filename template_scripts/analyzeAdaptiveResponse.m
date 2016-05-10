clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

dataDir = '/Volumes/Seagate Backup Plus Drive/160317_ibidi_RIactivin';
%dataDir = '/Users/idse/data_tmp/cycloheximide_after_20160330_32055 PM';
%dataDir = '/Users/idse/data_tmp/cycloheximide_before_20160330_42945 PM';
% LOOK INTO HOW THE LAST 10 TIME POINTS GOT LOST FOR CYCLOHEXIMIDE

meta = MetadataAndor(dataDir);
%meta.nTime = 110; % JUST FOR cycloheximide_after
%meta.nTime = 3; % JUST FOR cycloheximide_before
filenameFormat = meta.filename;

% manual metadata
%-------------------

% TODO : modify MetadataAndor to contain all info below

barefname = 'RIactivin100';
treatmentTime = 8; % first time point after treatment
conditions = {'RI + Activin 100 ng/ml'};
posPerCondition = 16;
nWells = 1;

% barefname = 'cycloheximide_after';
% %barefname = 'cycloheximide_before';
% treatmentTime = 4;
% conditions = {'no treatment','Activin 100 ng/ml', 'BMP 50 ng/ml', 'cyclohex 50 \mu g/ml',...
%               'MG','Activin + cyclohex','Activin + MG','BMP + cyclohex'};
% posPerCondition = 4;
% nWells = 8;

nucChannel = 2;
S4Channel = 1;

% visualize positions
%---------------------

% meta.displayPositions;

% TODO: create merged cellData for montage
% movies of distribution over time

%%

load(fullfile(dataDir,'positions'));

%csvwrite(fullfile(dataDir,'meanTimeTrace'),ratioMean);

%%

s = strsplit(meta.timeInterval,' ');
dt = str2double(s{1});
unit = s{2};
t = ((1:positions(1).nTime) - treatmentTime)*dt;

cd(dataDir);
saveResult = true;
minNCells = 10; % minimal number of cells
fgc = 'w';
bgc = 'k';

for wellnr = 1:nWells

    conditionPositions = posPerCondition*(wellnr-1)+1:posPerCondition*wellnr;

    ttraceCat = cat(1,positions.timeTraces);
    ttraceCat = ttraceCat(conditionPositions);
    nucMean = mean(cat(2,ttraceCat.nucLevelAvg),2);
    cytMean = mean(cat(2,ttraceCat.cytLevelAvg),2);
    bgMean = mean(cat(2,ttraceCat.background),2);

    clf 
    hold on
    ratioMean = (nucMean - bgMean)./(cytMean - bgMean);
    bad = any(cat(1,positions(conditionPositions).ncells) < minNCells,1);
    ratioMean(bad) = NaN;
    plot(t, ratioMean,'w','LineWidth',2)

    fs = 24;
    xlabel(['time (' unit ')'], 'FontSize',fs,'FontWeight','Bold','Color',fgc)
    ylabel('nuclear : cytoplasmic Smad4', 'FontSize',fs,'FontWeight','Bold','Color',fgc);

    axis([t(1), t(end)+50, 0.4, 1.5]);
    set(gcf,'color',bgc);
    set(gca, 'LineWidth', 2);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca,'XColor',fgc);
    set(gca,'YColor',fgc);
    set(gca,'Color',0.5*[1 1 1]);

    c = 1/80;
    plot(t,0.5+ (t>0).*t.*exp(-c*t)/50,'LineWidth',2)
    
    frame = export_fig(gcf,'-native -m2');
    
    hold off
end




