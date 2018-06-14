function [radialAvg, r] = makeAveragesNoSegmentation(...
                datadir, colSize, DAPIChannel, filenr, doubleNormalize)

% doubleNormalize: boolean
% first normalize by DAPI, then scale all profiles from 0 to 1

if ~exist('DAPIChannel','var')
    DAPIChannel = 1;
end
if ~exist('doubleNormalize','var')
    doubleNormalize = false;
end
if ~exist('filenr','var')
    filenr = [];
else
    filenr = ['_' num2str(filenr)];
end

load(fullfile(datadir,['colonies' filenr '.mat']));
load(fullfile(datadir,'metaData.mat'));

if ~exist('colSize','var')
    colSize = colonies(1).radiusMicron;
end

%restrict to colonies of the correct size
inds = [colonies.radiusMicron] == colSize;
colonies = colonies(inds);

r = imfilter(colonies(1).radialProfile.BinEdges,[1 1]/2)*meta.xres;
r(1) = colonies(1).radialProfile.BinEdges(1)*meta.xres;
r = r(1:end-1);
colCat = cat(3,colonies(:).radialProfile);

nucAvgAll = mean(cat(3,colCat.NucAvg),3);
nucAvgAllNormalized = bsxfun(@rdivide, nucAvgAll, nucAvgAll(:,DAPIChannel));

if doubleNormalize
    
    % make a version scaled from 0 to 1
    norm = max(nucAvgAllNormalized) - min(nucAvgAllNormalized);
    nucAvgDoubleNormalized = bsxfun(@minus, nucAvgAllNormalized, min(nucAvgAllNormalized));
    nucAvgDoubleNormalized = bsxfun(@rdivide, nucAvgDoubleNormalized', norm')';
    
    radialAvg = nucAvgDoubleNormalized;    
else
    radialAvg = nucAvgAllNormalized;
end


