%% locate code and data
addpath(genpath('C:\Users\mcg3\Documents\GitHub\stemcells'));
dataDir = 'F:\160627-QuickWashvsHourlyWash';
outputDir = 'F:\160627-QuickWashvsHourlyWash';
barefname = 'quickwash';
cd(dataDir);

numPos = 18;
numChannels = 2;

%% combine input files into output

for fi = 12:numPos
    for ci = 1:numChannels
        disp(['processing position ',num2str(fi),' channel ',num2str(ci-1)]);
        infiles = {...
            ['Run1\',sprintf([barefname '_f%.4d_w%.4d','.tif'],fi-1,ci-1)],...
            ['Run2\',sprintf([barefname '_f%.4d_t0000_w%.4d','.tif'],fi-1,ci-1)],...
            ['Run2\',sprintf([barefname '_f%.4d_t0001_w%.4d','.tif'],fi-1,ci-1)]};
            outfilebase = [sprintf([barefname '_f%.4d_w%.4d','.tif'],fi-1,ci-1)];
            outfile = fullfile(outputDir, outfilebase);
            catTiff(infiles,outfile);
    end
end