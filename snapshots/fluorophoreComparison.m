addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 

[scriptPath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;
cd(dataDir);

lasers = [405; 445; 488; 514; 561; 637];
filters = [450; 480; 525; 540; 600; 700];
          
Nlasers = numel(lasers);
Nfilters = numel(filters);

laseridx = 1:Nlasers;
filteridx = 1:Nfilters;

name = 'mTagRFP;Smad4GFP';
defaultIdx = 5; % the right laser/filter combination

%%

% determine file names
listing = dir(fullfile(dataDir,'*ims'));
fnames = {};
for fi = 1:Nfilters
    for li = 1:Nlasers
        substr = [num2str(lasers(li)) 'ex ' num2str(filters(fi)) 'em'];
        for k = 1:numel(listing)
            if contains(listing(k).name, substr)
                fnames{fi,li} = listing(k).name;
            end
        end
    end
end

% load images
imgs = {};
imgsMIP = {};
for fi = 1:Nfilters
    for li = 1:Nlasers
        if ~isempty(fnames{fi,li})
            [img, omeMeta] = readStack2(fnames{fi,li});
            imgs{fi,li} = img(:,:,:,3);
            imgsMIP{fi,li} = max(img,[],4);
        end
    end
end
%%
% determine limits
tol = 0.01;
minVals = nan([Nlasers Nfilters]);
maxVals = nan([Nlasers Nfilters]);
for fi = 1:Nfilters
    for li = 1:Nlasers
        if ~isempty(fnames{fi,li})
            lowhigh = stretchlim(imgsMIP{fi,li}, tol);
            minVals(fi,li) = lowhigh(1);%min(slice(:));
            maxVals(fi,li) = lowhigh(2);%max(slice(:));
        end
    end
end

minAll = min(minVals(:));
maxAll = max(maxVals(:));
tolAll = [minVals(defaultIdx,defaultIdx) maxVals(defaultIdx,defaultIdx)]
%%
% smaller previews to copy to slide
tolAll = [0 0.05];
margin = 25;
fs = 14;

screensize = get( 0, 'Screensize' );
w = screensize(3);
h = Nfilters*(screensize(3)/m + margin/2);
figure('Position', [1, 1, w, h]);

for fi = 1:Nfilters
    for li = 1:Nlasers

        i = (fi-1)*Nlasers + li;
        
        if ~isempty(imgs{fi,li})
            subplot_tight(Nfilters, Nlasers,i)
            imshow(imadjust(imgsMIP{fi,li}, tolAll));
            labelstr = ['\color{red}' 'ex' num2str(lasers(li)) '; em ' num2str(filters(fi))];
            text(margin, size(slice,1) - 2.5*margin, labelstr,'FontSize',fs,'FontWeight','bold');
        end
    end
end
fname = fullfile(dataDir,['crosstalk_' name '.png']);
saveas(gcf, fname);
%close;
