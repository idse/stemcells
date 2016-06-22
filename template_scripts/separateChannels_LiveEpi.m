%% locate code and data
addpath(genpath('C:\Users\Thomas\Documents\GitHub\stemcells')); 
cd('F:\Clayton\160620-ManualPulseTest')
%oifToTiffTimeLapse('Image966.vsi','D:\160602_pulsetest1\TIFF_FILES\I966');

%% separate files
for i = 20:35
    infile = ['Process_' num2str(i) '.vsi'];
    outfilebase = ['F:\Clayton\160620-ManualPulseTest\MIP\ManualPulseTest_MIP_p00' num2str(i-20) '_w0000'];
    oifToTiffTimeLapse(infile,outfilebase)
    disp(i)
end
