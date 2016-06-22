%% locate code and data
addpath(genpath('C:\Users\Thomas\Documents\GitHub\stemcells'));
dataDir = 'F:\Clayton\160621-ManualPulseTest2';
barefname = 'ManualPulseTest2';
cd(dataDir);

%% separate files
for i = 20:35
    infile = ['Process_' num2str(i) '.vsi'];
    outfilebase = fullfile(dataDir, sprintf([barefname '_MIP_p%.4d'],i-20));
    oifToTiffTimeLapse(infile,outfilebase)
    disp(i)
end
