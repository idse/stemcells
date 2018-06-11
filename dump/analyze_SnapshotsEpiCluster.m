clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

dataDir = '/Users/Idse/data_tmp/cellfate/pulses/170222_SyringePumpAendo2';
filelist = {'Process_1918.vsi','Process_1919.vsi','Process_1920.vsi','Process_1921.vsi'};

vsifile = fullfile(dataDir,filelist{1});
[~,barefname,~] = fileparts(vsifile);

resultsDir = fullfile(dataDir,'results');
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end

% TODO : local normalization by DAPI

% metadata
%----------------

meta = Metadata(vsifile);
    
% manually entered metadata
%------------------------------

meta.channelLabel = {'Sox17','H2B','Smad4','Sox2'};
meta.conditions = {'mtesr','+A10', '+A100', '+A10-1hr-pulses'};

nucChannel = 2;
dataChannels = [2 1 4];

markerChannels = setdiff(dataChannels, nucChannel);

% set manual intensity limits so images can be compared
Ilim = {[100 2^16-1],[100 1500],[100 2000],[100 500]};

%% load data for parameter check

fi = 1;
vsifile = fullfile(dataDir,filelist{fi});
[~,barefname,~] = fileparts(vsifile);

n = 1;
m = 2;
ymin = n*2048;
xmin = n*2048;
ymax = ymin + m*2048;
xmax = xmin + m*2048;
preview = {};
preview8bit = {};

series = 1;
img_bf = bfopen_mod(vsifile,xmin,ymin,xmax-xmin+1,ymax-ymin+1,series);

for ci = 1:numel(dataChannels)
    
    preview{ci} = img_bf{1}{dataChannels(ci),1};                    
    preview8bit{ci} = uint8((2^8-1)*mat2gray(preview{ci}, Ilim{dataChannels(ci)}));
end

previewChannels = 1:3;
filename = fullfile(dataDir, 'preview', [barefname '_RGBpreview.tif']);
imwrite(cat(3,preview8bit{previewChannels}), filename);
imshow(cat(3,preview8bit{previewChannels}))

%% load segmentation to test parameters

vsifile = fullfile(dataDir, filelist{fi});
P = Position(meta.nChannels, vsifile, meta.nTime);
P.setID(pi);
seg = P.loadSegmentation(fullfile(dataDir,'MIP'), nucChannel);
segpart = seg(xmin:xmax, ymin:ymax);

%% finding the right options for extract data

% externally I will have all indices starting at 1
% the Andor offset to start at 0 will be internal

opts = struct(  'cytoplasmicLevels',    false,... 
                    'dataChannels',     1:meta.nChannels,... 
                    'tMax',             meta.nTime,...
                    'segmentationDir',  fullfile(dataDir,'MIP'),...
                    'nucShrinkage',     2,...
                    'cytoSize',         5,...
                    'cytoMargin',       5,...
                    'bgMargin',         4);
                
opts.cleanupOptions = struct('separateFused', true,'openSize',5,...
    'clearBorder',true, 'minAreaStd', 1, 'minSolidity',0, 'minArea',500);

% try out the nuclear cleanup settings on some frame:
segpartclean = nuclearCleanup(segpart, opts.cleanupOptions);

% visualize
s = 0.2;
nucim = mat2gray(preview8bit{dataChannels==nucChannel});

RGBsegCheck = cat(3,   nucim + s*mat2gray(segpartclean),...
                        nucim,...
                        nucim + s*segpart);
figure, imshow(RGBsegCheck)

filename = fullfile(resultsDir, [barefname 'segCheck.tif']);
imwrite(RGBsegCheck, filename);

%% process data for first image

debugInfo = P.extractData(dataDir, nucChannel, opts);

nucmaskpart = debugInfo.nucmask(xmin:xmax, ymin:ymax);
figure, imshow(cat(3,   nucim + s*nucmaskpart,...
                        nucim,...
                        nucim + s*segpart))
                    
% save
save(fullfile(resultsDir, [barefname '_positions.mat']),'P');

%% scatter cell locations on top of part of preview to check result

XY = P.cellData.XY;

inpreview =     XY(:,1) < xmax & XY(:,1) > xmin...
                & XY(:,2) < ymax & XY(:,2) > ymin;
XY = XY(inpreview,:);

imshow(nucim,[])
hold on
scatter(XY(:,1) - xmin,XY(:,2) - ymin, 200, '.r');
hold off

%% process the other images

for fi = 2:numel(filelist)
    
    vsifile = fullfile(dataDir, filelist{fi});
    [~,barefname,~] = fileparts(vsifile);

    P = Position(meta.nChannels, vsifile, meta.nTime);
    P.setID(pi);
    P.extractData(dataDir, nucChannel, opts);

    save(fullfile(resultsDir, [barefname '_positions.mat']),'P');
end

%% make distributions

% make distributions out of processed data
stats = getStats(dataDir, filelist); % not normalized
% stats = getStats(dataDir, filelist, nucChannel); % normalized
stats.channelLabel = meta.channelLabel;
stats.conditions = meta.conditions;

%% plot distributions

% overlay of the distribution of different conditions
for channelIndex = markerChannels

    plotDistribution(stats, channelIndex)

    filename = fullfile(resultsDir,...
            ['distOverlay_' meta.channelLabel{channelIndex}]);
    saveas(gcf, filename);
    saveas(gcf, [filename '.png']);
    
    cumulative = true;
    plotDistribution(stats, channelIndex, cumulative)

    filename = fullfile(resultsDir,...
            ['cumdistOverlay_' meta.channelLabel{channelIndex}]);
    saveas(gcf, filename);
    saveas(gcf, [filename '.png']);
end

%% scatter plots of markers against each other

%markerChannels = setdiff(dataChannels, nucChannel);
fs = 15;

for fi = 1:numel(filelist)
    
    for ci1 = 1:numel(markerChannels);
        for ci2 = ci1+1:numel(markerChannels)

            c1 = markerChannels(ci1);
            c2 = markerChannels(ci2);
            
            A = stats.nucLevel{fi}(:,c1);
            B = stats.nucLevel{fi}(:,c2);

            scatter(A, B,1,'.');
            xlabel(meta.channelLabel{c1}, 'FontSize',fs, 'FontWeight','Bold');
            ylabel(meta.channelLabel{c2}, 'FontSize',fs, 'FontWeight','Bold');
            title([meta.conditions{fi} '      (' filelist{fi} ')'],...
                'Interpreter','none', 'FontSize',fs, 'FontWeight','Bold');
            
            axis([stats.lim{c1}' stats.lim{c2}']);

            set(gcf,'color','w');
            set(gca, 'LineWidth', 2);
            set(gca,'FontSize', fs)
            set(gca,'FontWeight', 'bold')

            filename = fullfile(resultsDir, ['scatter_' meta.channelLabel{c1} '_' meta.channelLabel{c2} '_' barefname]);
            saveas(gcf, filename);
            saveas(gcf, [filename '.png']);
        end
    end
end


%% try clustering in 2D

k = 2;
fi = 1;
c1 = 1; c2 = 4;
X = [stats.nucLevel{fi}(:,c1) stats.nucLevel{fi}(:,c2)];
for fi = 2:4
    X = cat(1,X,[stats.nucLevel{fi}(:,c1) stats.nucLevel{fi}(:,c2)]);
end
X(X(:,2) > stats.lim{c1}(1),:) = [];
%scatter(X(:,1), X(:,2), 1, '.');

Xnorm = X;
for i = 1:2
    Xnorm(:,i) = Xnorm(:,i)./max(Xnorm(:,i));
end

% % % kmeans -> doesn't work well
%y = kmeans(X,k);

% Gaussian mixture

gmfit = fitgmdist(Xnorm,k,'CovarianceType','full','SharedCovariance',false,...
        'RegularizationValue', 10^(-3));
%gmfit = fitgmdist(X,k,'CovarianceType','full','SharedCovariance',false);
[y, ~, P, logpdf] = cluster(gmfit,Xnorm);

%idx = logpdf>-10;
cutoff = 0.95;
idx = P(:,1) > cutoff | P(:,2) > cutoff;
scatter(Xnorm(:,1), Xnorm(:,2), 1,'.');

scatter(Xnorm(idx,1), Xnorm(idx,2), 1, y(idx),'.');

hold on
% means of gaussians
scatter(gmfit.mu(:,1), gmfit.mu(:,2),100,'.r')

% covariance
for i = 1:2
    
    [V,D] = eig(gmfit.Sigma(:,:,i));
    t = linspace(0,2*pi);
    
    % major and minor axis length
    lmin = sqrt(D(1,1));
    lmaj = sqrt(D(2,2));

    majV = V(:,2);
    phi = atan(majV(2)/majV(1));
    R = [cos(phi) sin(phi); -sin(phi) cos(phi)];

    %quiver(gmfit.mu(i,1), gmfit.mu(i,2), majV(1), majV(2));
    
    X = lmaj*cos(t);
    Y = lmin*sin(t);
    
    Xp = R(1,1)*X + R(2,1)*Y + gmfit.mu(i,1);
    Yp = R(1,2)*X + R(2,2)*Y + gmfit.mu(i,2);
    
    plot(Xp, Yp,'LineWidth',2)
end
hold off

%% per condition

fi = 2;
Xc = [stats.nucLevel{fi}(:,c1) stats.nucLevel{fi}(:,c2)];
size(Xc)
Xcnorm = Xc;
for i = 1:2
    Xcnorm(:,i) = Xcnorm(:,i)./max(Xcnorm(:,i));
end
[y, nlogl, P, logpdf] = cluster(gmfit,Xcnorm);
% idx = logpdf>-10;
colormap lines(2);

idx = P(:,1) > cutoff | P(:,2) > cutoff;
scatter(Xc(idx,1), Xc(idx,2), 1, y(idx),'.');
xlim(stats.lim{c1}');
ylim(stats.lim{c2}');

%% try hierarchical clustering -> scales poorly

N = 10000;
idx = round((length(X)-1)*rand([N 1])+1);
Xred = X(idx,:);
for i = 1:2
    Xred(:,i) = Xred(:,i)./max(Xred(:,i));
end

% Z = linkage(Xred,'ward','euclidean');
% T = cluster(Z,'maxclust',2);

Z = linkage(Xred,'average','euclidean');
T = cluster(Z,'cutoff',0.4,'criterion','distance');%
scatter(Xred(:,1),Xred(:,2),20,T,'.')

% T = cluster(Z,'maxclust',2);
bins = unique(T);
n = hist(T,bins);

outliers = bins(n<20);
for i = 1:numel(outliers)    
    Xred(T == outliers(i),:) = [];
    T(T == outliers(i),:) = [];
end

%dendrogram(Z)
clf
scatter(Xred(:,1),Xred(:,2),20,T,'.')

%% conditional prob based on cluster

fi = 2;
ci = 2;
X = [stats.nucLevel{fi}(:,c1) stats.nucLevel{fi}(:,c2)];
y = mixGaussPred(X',model);
% scatter(X(:,1), X(:,2), 1, y,'.');
% axis([stats.lim{c1}' stats.lim{c2}']);

bins = linspace(stats.lim{c2}(1),stats.lim{c2}(2),50);
[n,bins] = hist(X(:,ci),bins);
[n1,bins] = hist(X(y==1,ci),bins);
n2 = hist(X(y==2,ci),bins);
n1 = n1./(sum(n1)+sum(n2));
n2 = n2./(sum(n1)+sum(n2));
n = n./sum(n);
plot(bins, n)
hold on
plot(bins, n1,'g')
plot(bins, n2,'r')
hold off

%% conditional probabilities

fi = 1;
ci = 1;
n = stats.histograms{fi, ci};
bins = stats.bins{fi,ci};
dist = n./sum(n);
cumdist = cumsum(n./sum(n));
i = find(cumdist > 0.99, 1, 'first');
thresh = bins(i);

%%
fi = 1;
ci = 2;

X = [stats.nucLevel{fi}(:,c1) stats.nucLevel{fi}(:,c2)];
y = 1 + (X(:,1) > thresh);

bins = linspace(stats.lim{c2}(1),stats.lim{c2}(2),50);
[n,bins] = hist(X(:,ci),bins);

[n1,bins] = hist(X(y==1,ci),bins);
n2 = hist(X(y==2,ci),bins);
n1 = n1./(sum(n1)+sum(n2));
n2 = n2./(sum(n1)+sum(n2));
n = n./sum(n);
plot(bins, n)
hold on
plot(bins, n1,'g')
plot(bins, n2,'r')
hold off

%% determine fraction of sox17+ etc

cutoff = 2;
channelIdx = 2;

for fi = 1:numel(filelist)
    
    p = allData{fi};
    
    N = p.cellData.nucLevel(:,nucChannel);
    A = p.cellData.nucLevel(:,channelIdx)./N;

    disp([filelist{fi} ' : ' num2str(100*sum(A > cutoff)/numel(A),2)]);
end
