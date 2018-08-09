function output = makeAveragesNoSegmentation(...
                meta, colSize, DAPIChannel, colonies)

% doubleNormalize: boolean
% first normalize by DAPI, then scale all profiles from 0 to 1

if ~exist('DAPIChannel','var')
    DAPIChannel = 1;
end

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
if ~isempty(DAPIChannel)
    nucAvgAllDAPINormalized = bsxfun(@rdivide, nucAvgAll, nucAvgAll(:,DAPIChannel));

    % make a version scaled from 0 to 1
    norm = max(nucAvgAllDAPINormalized) - min(nucAvgAllDAPINormalized);
    nucAvgDoubleNormalized = bsxfun(@minus, nucAvgAllDAPINormalized, min(nucAvgAllDAPINormalized));
    nucAvgDoubleNormalized = bsxfun(@rdivide, nucAvgDoubleNormalized', norm')';
else
    nucAvgAllDAPINormalized = [];
    nucAvgDoubleNormalized = [];
end

output = struct('nucAvg', nucAvgAll,...
                'nucAvgDAPINormalized', nucAvgAllDAPINormalized,...
                'nucAvgDAPImaxNormalized', nucAvgDoubleNormalized,...
                'r',r,...
                'colSize',colSize);
            

