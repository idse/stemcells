function avgs = nucAvgsNoSegment(dataDir)
%return nuclear average from colonies.mat

load(fullfile(dataDir,'colonies.mat'));
load(fullfile(dataDir,'metaData.mat'));
DAPIChannel = find(strcmp(meta.channelNames,'DAPI'));

if isempty(DAPIChannel)
DAPIChannel = 1;
end


coloniesToDo = colonies;
colCat = cat(3,coloniesToDo(:).radialProfile);

nucAvgAll = mean(cat(3,colCat.NucAvg),3);
%cytAvgAll = mean(cat(3,colCat.CytAvg),3);

nucAvgAllNormalized = bsxfun(@rdivide, nucAvgAll, nucAvgAll(:,DAPIChannel));

avgs = nucAvgAllNormalized;