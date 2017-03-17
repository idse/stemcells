clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/')); 
%dataDir = '/Users/idse/data_tmp/161107_NODALcrispr';
dataDir = '/Volumes/IdseData/170210_ANBMP4withSB/part1';

nucChannel = 1;
S4Channel = 0;

channels = [nucChannel S4Channel];

% save idx for channels containing nuclear marker if multiple z-slices 
saveidx = [true false]; 

inputdir = dataDir;
outputdir = fullfile(dataDir, 'MIP');

%% 

% laser
%---------

if exist(fullfile(dataDir, 'MATL_Mosaic.log'),'file')
    
    batchMIP_lsm(inputdir, outputdir, channels, saveidx);
end

% epi 
%---------
if ~isempty(dir(fullfile(dataDir, '*.vsi')))
    
    batchMIP_epi(inputdir, outputdir, channels, saveidx); %, meta.nTime, 1);
    
% Andor 
%---------
elseif ~isempty(dir(fullfile(dataDir, '*.txt')))
    
    outputdir = fullfile(dataDir, 'MIP');
    batchMIP_Andor(inputdir, outputdir, channels, saveidx);
%     
%     meta = MetadataAndor(dataDir);
%     type = 'SIP';
%     batchMIP_Andor(inputdir, outputdir, channels, saveidx,...
%         meta.nTime, 1:meta.nPositions, type);
end
