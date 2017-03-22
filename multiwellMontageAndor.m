clear all;
addpath(genpath('/Users/idse/repos/Warmflash/stemcells'));

% dataDir = '/Users/idse/data_tmp/';
% filename = fullfile(dataDir, 'IH_96wellTest.xyz');

dataDir = '/Volumes/STORENGO';
filename = fullfile(dataDir,'160317_IH.xyz');

resolution = 0.325; % 40x
micropattern  = false;


%% read positions from XYZ file

positions = struct();
fid = fopen(filename);

if fid == -1
    error('file not found');
end

% check first line
tline = fgets(fid);
if ~strcmp(strtrim(tline),'[XYZScan]')
    error('not an XYZ file');
end

% scan through the header
tline = fgets(fid);
while ischar(tline) 
    k = strfind(tline,'=');
    if isempty(k)
        break;
    else
        val = strtrim(tline(k + 1:end));
        fieldname = strrep(tline(1:k-1),' ', '_');
        positions.(fieldname) = val;
    end
    tline = fgets(fid);
end
fclose(fid);

% extract positions properly
s = strsplit(positions.XYFields,'\t');
nPos = str2double(s{1});
XYZ = zeros([nPos 3]);

for i = 1:nPos
    strsplit(s{i+1},',');
    XYZ(i,:) = str2double(strsplit(s{i+1},','));
end
%positions.XYFields = XYZ;

if ~isempty(positions.AutofocusPositions)
    s = strsplit(positions.AutofocusPositions,'\t');
    autofocus = zeros([nPos 1]);

    for i = 1:nPos
        strsplit(s{i+1},',');
        autofocus(i,:) = str2double(s{i+1});
    end
end
%positions.AutofocusPositions = autofocus;

%% make grid out of single position (or average of two)

% 3x3 grid will have length (2*(1024-overlap) + 1024-2*overlap)*resolution
% at 40x that is 868 micron so large enough for a 700 micron colony
overlapPixel = 100; 

gridSize = [5 5]; % ODD NUMBERS HERE FOR NOW
nGridPositions = gridSize(1)*gridSize(2);

gridXYZ = {};
AF = {};
if micropattern
    nColonies = size(XYZ,1)/2;
else
    nColonies = size(XYZ,1);
end

for pi = 1:nColonies

    if micropattern
        x = mean(XYZ(2*pi-1:2*pi,:),1);
        af = mean(autofocus(2*pi-1:2*pi));
    else
        x = XYZ(pi,:);
        af = autofocus(pi);
    end
    
    spacing = (1024 - overlapPixel)*resolution;

    gridXYZ{pi} = zeros([nGridPositions 3]);
    AF{pi} = zeros([nGridPositions 1]);

    nmax = (gridSize(1)-1)/2;
    mmax = (gridSize(2)-1)/2;
    
    i = 1;
    for n = -nmax:nmax
        for m = -mmax:mmax
            gridXYZ{pi}(i,:) = [x(1) + m*spacing, x(2) + n*spacing, x(3)];
            AF{pi}(i) = af;
            i = i+1;
        end
    end
end
%scatter(gridXYZ{1}(:,1), gridXYZ{1}(:,2));

% convert XYZ back to string
gridXYZcombined = cat(1,gridXYZ{:});
AFCombined = cat(1,AF{:});
newNPos = size(gridXYZcombined,1);
posStrings = {};
AFstrings = {};
for i = 1:newNPos
    AFstrings{i} = num2str(AFCombined(i));
    posStrings{i} = sprintf('%d,%d,%.2f',gridXYZcombined(i,:));
end
newXYFields = strjoin([num2str(newNPos) posStrings],'\t');
newAF = strjoin([num2str(newNPos) AFstrings],'\t');

newPositions = positions;
newPositions.XYFields = newXYFields;
newPositions.AutofocusPositions = newAF;

%% check

scatter(gridXYZcombined(:,1), gridXYZcombined(:,2),'g');
hold on
scatter(XYZ(:,1),XYZ(:,2),'r');
hold off
axis equal

%% write back to file

[~,barefname,~] = fileparts(filename);
outFilename = fullfile(dataDir, [barefname '_montage.xyz']);
fileID = fopen(outFilename,'w');

fields = fieldnames(newPositions);

fprintf(fileID,'[XYZScan]\n'); 

for i = 1:numel(fields)
    
    line = [fields{i} '=' newPositions.(fields{i}),'\n'];
    fprintf(fileID,line); 
end

fclose(fileID);

