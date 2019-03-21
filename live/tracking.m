clear all; close all;

addpath(genpath('/Users/idseimac/stemcells')); 

%[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);

dataDir = '/data/170705_Gridse/MIP';

%%
ti = 1;
fname = sprintf('%0.5d.h5',ti);
h5disp(fullfile(dataDir, 'Gridse_MIP_p0000_w0001_H5-Event-Sequence', fname));

%%
data = h5read(fullfile(dataDir, 'Gridse_MIP_p0000_w0001_H5-Event-Sequence', fname), '/tracking/Moves');

%%

h5disp(fullfile(dataDir,'Gridse_MIP_p0000_w0001_Simple Segmentation.h5'))