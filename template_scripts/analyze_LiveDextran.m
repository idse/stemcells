clear all; close all;

% directories
addpath(genpath('C:\Users\Thomas\Documents\GitHub\stemcells'));
dataDir = 'D:\160630_fluorDextran40xBinning';

% metadata
meta = MetadataAndor(dataDir);
filenameFormat = meta.filename;

% other information
barefname = 'fluorDextran40xBinning';
posPerCondition = 10;
nWells = 8;
bins = 16;
conditions = {'00','10','20','30','40','50','60','70'};
units = 'ng/mL';
concentrations = [0,10,20,30,40,50,60,70];

% manually entered meta data
meta.xSize = meta.xSize/sqrt(bins);
meta.ySize = meta.ySize/sqrt(bins);

% visualize positions
%meta.displayPositions;

%% extract data

flipOrder = false; %flip the order of the data for the last 4 wells

% dimensions are ypixel, xpixel, well position, well, image cycle
positions = zeros(meta.ySize,meta.xSize,posPerCondition,nWells,meta.tPerFile);

% create positions data file
for ti = 1:meta.tPerFile
    for wi = 1:nWells
        for pi = 1:posPerCondition
            timePoint = (ti-1)*meta.nPositions + (wi-1)*posPerCondition + pi;
            if mod(timePoint,nWells*posPerCondition) == 0
                fprintf('\n');
            else
                fprintf('.');
            end
            positions(:,:,pi,wi,ti) = ...
                imread(fullfile(dataDir,filenameFormat),timePoint);
        end
    end
end

% correct for different imaging order
if nWells > 4 && flipOrder
    positions(:,:,5:end,:,:) = positions(:,:,end:-1:5,:,:);
end

disp('Saving data...');
save(fullfile(dataDir,'positions'), 'positions');
disp('Done');

%% center subset of the data

index = [(meta.xSize*1/4):(meta.xSize*3/4);(meta.ySize*1/4):(meta.ySize*3/4)];
positions = positions(index(2,:),index(1,:),:,:,:);

%% spatial background subtraction

% find min spatial values for background
b = squeeze(min(min(min(positions(:,:,:,:,:),[],3),[],4),[],5));

% subtract background spatials
s = size(positions); s(1:2) = [1,1];
positions = positions - repmat(b,s);

%saveTiff(positions,fullfile(dataDir,'bgSubtraction.tif'));

%% segment out bright spots

save_tiff = true;  %save a tiff file of the segmentation
throw_out = true;  %throw out images with large bright spots

seg_diskSize = 3;  %large disk sizes will make processing slow
seg_stdDev = 2;    %large stdDev may not capture all bright spots
seg_minPix = 8;    %min pixel size of segments to throw out small spots
throw_size = 120;  %min pixel size of segmentation to throw out images

% segmentation
for ti = 1:meta.tPerFile
    for wi = 1:nWells
        for pi = 1:posPerCondition
            timePoint = (ti-1)*meta.nPositions + (wi-1)*posPerCondition + pi;
            % segment and index
            img = squeeze(positions(:,:,pi,wi,ti));
            seg = imseg2(img,seg_diskSize,seg_stdDev,seg_minPix);
            if sum(seg(:)) > throw_size %position is compromised
                positions(:,:,pi,wi,ti) = NaN;
                throw = true;
            else
                img(seg) = NaN;
                positions(:,:,pi,wi,ti) = img;
                throw = false;
            end
            % progress indicator
            if mod(timePoint,posPerCondition) == 1
                fprintf(['Time ' num2str(ti) ', Well ' num2str(wi) ': ']);
            end
            if throw
                fprintf('/');
            else
                fprintf('.');
            end
            if mod(timePoint,posPerCondition) == 0
                fprintf('\n');
            end
            % save segmentation into tiff file
            if save_tiff
                if ~any([ti,wi,pi]-1) %first image
                    imwrite(seg,fullfile(dataDir,'segmentedSpots.tif'),'Compression','none');
                else
                    imwrite(seg,fullfile(dataDir,'segmentedSpots.tif'),'writemode','append','Compression','none');
                end
            end
        end
    end
end

saveTiff(positions,fullfile(dataDir,'segmentedData.tif'));

%% spatial rescaling
% make sure this step is after bg sub, segmentation, etc.

% find the rescaling matrix
pixelMeans = squeeze(mean(mean(mean(positions(:,:,:,:,:),3,'omitnan'),4,'omitnan'),5,'omitnan'));
xScale = mean(pixelMeans,1); xScale = xScale / max(xScale);
yScale = mean(pixelMeans,2); yScale = yScale / max(yScale);
r = yScale * xScale;
imshow(r); %visualize the scaling

% rescale positions
s = size(positions); s(1:2) = [1,1];
positions = positions ./ repmat(r,s);

saveTiff(positions,fullfile(dataDir,'rescaledData.tif'));

%% spatial rescaling by well

s = size(positions); s([1,2,4]) = [1,1,1];
for wi = 1:nWells
    pixelMeans = squeeze(mean(mean(positions(:,:,:,wi,:),3,'omitnan'),5,'omitnan'));
    xScale = mean(pixelMeans,1); xScale = xScale / max(xScale);
    yScale = mean(pixelMeans,2); yScale = yScale / max(yScale);
    r = yScale * xScale;
    %figure, imshow(r);
    positions(:,:,:,wi,:) = positions(:,:,:,wi,:) ./ repmat(r,s);
end

saveTiff(positions,fullfile(dataDir,'rescaledData.tif'));

%% flatten the pixels of each image  
    
pos = size(positions);
flatPix = reshape(positions,[pos(1)*pos(2),pos(3),pos(4),pos(5)]);
flatPix(:,:,[1,2,3,5,8,7,6,4],:) = flatPix(:,:,1:8,:);

disp('Saving data...');
save(fullfile(dataDir,'flatPix'), 'flatPix');
disp('Done');

%% load data

disp('Loading file...');
load(fullfile(dataDir,'positions'));
disp('Done');

disp('Loading file...');
load(fullfile(dataDir,'flatPix'));
disp('Done');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STATISTICAL ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find statistical results
% calculates these statistics for each position, then takes average by well.
% in general, I first evaluate each image, then take the median of the
% positions, and then take the mean across time points.

numStats = 3; %the number of statistics being calculated. Change this
              %if you add another evaluation.
wellsWanted = 1:8; %save time by processing only some wells
avgT = false; %average the statistics across time

% manually ignore positions. Row = well, Col = positions.
% we ignore across time points, because problems associated with one
% positions seem to persist across time points.
ignore = [...
    0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0;];
[w,p] = find(ignore == 1);
tmp = flatPix;
for i = 1:numel(w)
    tmp(:,p(i),w(i),:) = NaN;
end
tmp = tmp(:,:,wellsWanted,:);

% always same size. Wells not evaluated will just have zeros.
fluorStats = zeros(nWells,meta.tPerFile,numStats); %holds the evaluations

% find the mean
fluorStats(wellsWanted,:,1) = squeeze(mean(mean(tmp,1,'omitnan'),2,'omitnan'));
% find the variance
fluorStats(wellsWanted,:,2) = squeeze(mean(var(tmp,0,1,'omitnan'),2,'omitnan'));
% find the mode
edges = 1:1:squeeze(max(tmp(:)));
bins = discretize(tmp,edges);
fluorStats(wellsWanted,:,3) = squeeze(mean(mode(bins,1),2,'omitnan'));
clear tmp

% average across time
if avgT
    fluorStats = squeeze(mean(fluorStats,2));
end

save(fullfile(dataDir,'fluorStats'), 'fluorStats');

%% load statistical results

load(fullfile(dataDir,'fluorStats'));

%% plot statistical results

wellsWanted = 1:8; %compare just these wells

% means
figure
gscatter(concentrations(wellsWanted)',fluorStats(wellsWanted,1),...
    conditions(wellsWanted)',[],[],25);
title('Concentration vs Mean Intensity');
% variance
figure
gscatter(concentrations(wellsWanted)',fluorStats(wellsWanted,2),...
    conditions(wellsWanted)',[],[],25);
title('Concentration vs Variance of Intensity');
% mode
figure
gscatter(concentrations(wellsWanted)',fluorStats(wellsWanted,3),...
    conditions(wellsWanted)',[],[],25);
title('Concentration vs Mode Intensity');

%% correct for bg signal and noise

tmp = fluorStats;
bgWell = 1; %the well to use for calculating background intensity

% true signal = raw signal - background signal - noise
% noise = sqrt(raw signal)

% find bg value
bgMean = tmp(bgWell,1) - sqrt(tmp(bgWell,1));
bgVar = tmp(bgWell,2) - sqrt(tmp(bgWell,2));

% find noise values
shotNoiseMeans = sqrt(tmp(:,1));
shotNoiseVars = sqrt(tmp(:,2));

% calculate the ratios
trueMeans = tmp(:,1) - bgMean - shotNoiseMeans;
trueVars = tmp(:,2) - bgVar - shotNoiseVars; % shot noise adds no variability?
trueRatios = trueMeans./trueVars;

% visualize
plot(trueRatios(2:8))
corr(trueRatios(2:8),(1:7)')
% positive trend ratios = means sub too big or vars sub too small
% negative trend ratios = means sub too small or vars sub too big

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DENSITY DISTRIBUTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% histogram density of each condition

wellsWanted = 1:8;
posWanted = 1:10;
tpWanted = 1:3;

percentRange = [0.1,99.9]; %percentage range of data to plot
ymax = (numel(tpWanted) * numel(posWanted) * meta.xSize * meta.ySize) / 70;

% set up the subplot grid
if numel(wellsWanted) == 1
    gridSize = [1,1];
else
    gridSize = [ceil(numel(wellsWanted)/2),2];
end

% plot the histograms
figure
for wi = 1:numel(wellsWanted)
    pixels = flatPix(:,posWanted,wellsWanted(wi),tpWanted);
    pixels = pixels(~isnan(pixels));
    % finding the outliers
    indexRange = round(numel(pixels)*percentRange/100);
    tmp = sort(pixels);
    fluorRange = [tmp(indexRange(1)),tmp(indexRange(2))];
    % plotting
    subplot(gridSize(1),gridSize(2),wi);
    hist(pixels,fluorRange(1):fluorRange(2))
    axis([fluorRange(1) fluorRange(2) 0 ymax])
    %histogram(wellPixels,'BinLimits',fluorRange);
    title([conditions{wellsWanted(wi)} [' ' units]]);
end

%% overlaid histogram densities of each condition

wellsWanted = 1:8;
posWanted = 1:10;
tpWanted = 1:1;

binStep = 1; %higher step = lower resolution
Imin = 0;    %min fluorescence
Imax = 400;  %max fluorescence
ymax = (numel(tpWanted) * numel(posWanted) * meta.xSize * meta.ySize) / 70;

% flatten across positions
if ~exist('multiWellPixels','var')
    multiWellPixels = shiftdim(flatPix(:,posWanted,:,tpWanted),2);
    multiWellPixels = multiWellPixels(:,:);
end

% plotting
hist(multiWellPixels(wellsWanted,:)',Imin:binStep:Imax);
axis([Imin Imax 0 ymax]);
legend(strcat(conditions(wellsWanted), [' ' units]));
title('PDF comparison of Dextran conditions')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CUMULATIVE DISTRIBUTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cumulative distributions of conditions

wellsWanted = 1:8;
posWanted = 1:10;
tpWanted = 1:1;

binStep = 1; %higher step = lower resolution
Imin = 0;    %min fluorescence
Imax = 400;  %max fluorescence

colors = hsv(numel(wellsWanted));

% plot the distributions
figure
hold on
for wi = 1:numel(wellsWanted)
    pixels = flatPix(:,posWanted,wellsWanted(wi),tpWanted);
    pixels = pixels(~isnan(pixels));
    n = histc(pixels,Imin:binStep:Imax);
    plot(cumsum(n)./sum(n),'LineWidth',2,'Color',colors(wellsWanted(wi),:));
end
hold off
xlim([0,size(n,1)]);
legend(strcat(conditions(wellsWanted), [' ' units]));
title('CDF comparison of Dextran conditions')

%% Cumulative distributions of positions

wellsWanted = 6:7;
posWanted = 1:10;
tpWanted = 1:3;

binStep = 1; %higher step = lower resolution
Imin = 0;    %min fluorescence
Imax = 400;  %max fluorescence

colors = hsv(numel(posWanted));
linetype = {'-',':','-.','--','-',':','-.','--','-',':'};

% plot the distributions
for wi = 1:numel(wellsWanted)
    figure
    hold on
    for pi = 1:numel(posWanted)
        pixels = flatPix(:,posWanted(pi),wellsWanted(wi),tpWanted);
        pixels = pixels(~isnan(pixels));
        n = histc(pixels,Imin:binStep:Imax);
        plot(cumsum(n)./sum(n),linetype{pi},'LineWidth',2,'Color',colors(posWanted(pi),:));
    end
    hold off
    xlim([0,size(n,1)]);
    legend(strread(num2str(posWanted),'%s'));
    title(['CDF comparison of positions for ',[conditions{wellsWanted(wi)},' ',units]])
end

%% Cumulative distributions of times

wellsWanted = 1:8;
posWanted = 1:10;
tpWanted = 1:3;

binStep = 1; %higher step = lower resolution
Imin = 0;    %min fluorescence
Imax = 400;  %max fluorescence

colors = hsv(numel(tpWanted));

% plot the distributions
for wi = 1:numel(wellsWanted)
    figure
    hold on
    for ti = 1:numel(tpWanted)
        pixels = flatPix(:,posWanted,wellsWanted(wi),tpWanted(ti));
        pixels = pixels(~isnan(pixels));
        n = histc(pixels,Imin:binStep:Imax);
        plot(cumsum(n)./sum(n),'LineWidth',2,'Color',colors(tpWanted(ti),:));
    end
    hold off
    xlim([0,size(n,1)]);
    legend(strread(num2str(tpWanted),'%s'));
    title(['CDF comparison of times for ',[conditions{wellsWanted(wi)},' ',units]])
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% KOLMOGOROV-SMIRNOV TESTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KS tests for time

wellsWanted = 1:8;
testPairs = combnk(1:meta.nTime,2);
testResults = zeros(nWells,size(testPairs,1));
numPixels = 100000;   % how many pixels to sample for each condition
replace = false;      % sample w/ or w/o replacement. May change result

for tpi = 1:size(testPairs,1)
        disp([...
            'Comparing time ',num2str(testPairs(tpi,1)),...
            ' and time ',num2str(testPairs(tpi,2))]);
        for wi = 1:numel(wellsWanted)
            wellidx = wellsWanted(wi);
            [~,p] = kstest2(...
                randsample(squeeze(flatPix(testPairs(tpi,1),wellidx,:)),numPixels),...
                randsample(squeeze(flatPix(testPairs(tpi,2),wellidx,:)),numPixels));
            testResults(wellidx,tpi) = p;
        end
end

disp('Test pairs:'); disp(testPairs);
disp('Test results:'); disp(testResults);