%% locate code and data
addpath(genpath('C:\Users\mcg3\Documents\GitHub\stemcells'));
dataDir = 'G:\Clayton\160621-ManualPulseTest2';
barefname = 'ManualPulseTest2';
cd(dataDir);

%% separate files
for i = 36%:43
    infile = ['Process_' num2str(i) '.vsi'];
    outfilebase = fullfile([dataDir '\MIP'], sprintf([barefname '_MIP_p%.4d'],i-36));
    oifToTiffTimeLapse(infile,outfilebase);
    disp(i)
end