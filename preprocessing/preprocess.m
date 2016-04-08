clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/')); 

dataDir = '/Users/idse/data_tmp/cycloheximide_after_20160330_32055 PM';
dataDir = '/Users/idse/data_tmp/cycloheximide_before_20160330_42945 PM';

nucChannel = 1;
S4Channel = 0;

[meta, rawmeta] = readMeta_Andor(dataDir);

%%
% visualize positions

displayPositions_Andor(meta);

%%

% ASSUME: files are split by position & wavelength

channels = [nucChannel S4Channel];
saveidx = [true false];
inputdir = dataDir;
outputdir = fullfile(dataDir, 'MIP');

batchMIP_Andor(inputdir, outputdir, channels, saveidx);
