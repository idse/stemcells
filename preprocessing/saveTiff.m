function saveTiff(data,outfile)
% data must be array of five dimensions, with last three dims being
% position, wells, and time...in that order.

disp('Saving tiff...');

% find the light value range
cleanData = sort(data(~isnan(data)));
numPixels = numel(cleanData);
Imin = cleanData(round(0.02 * numPixels));
Imax = cleanData(round(0.98 * numPixels));

% write to file
S = size(data);
for ti = 1:S(5)
    for wi = 1:S(4)
        for pi = 1:S(3)
            if all([ti,wi,pi] == 1) %first image
                imwrite(mat2gray(data(:,:,pi,wi,ti),[Imin Imax]),...
                    outfile,'Compression','none');
            else
                imwrite(mat2gray(data(:,:,pi,wi,ti),[Imin Imax]),...
                    outfile,'Compression','none','writemode','append');
            end
        end
    end
end

disp('Done');

return
