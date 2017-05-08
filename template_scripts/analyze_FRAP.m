clear all; close all;

clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/stemcells')); 
dataDir = '/Users/idse/data_tmp/0_kinetics/170302_FRAPagain';

%oibfile = 'untreated.oib';
%oibfile = 'LMB1h.oib';
%oibfile = 'LMB4h.oib';
%oibfile = 'ALMB2h20min.oib';
oibfile = 'notreatNextDay.oib';

[~,barefname,~] = fileparts(oibfile);

nucChannel = 2;
S4Channel = 1;

r = bfGetReader(fullfile(dataDir,oibfile));
omeMeta = r.getMetadataStore();

tres = double(omeMeta.getPixelsTimeIncrement(0).value);

%% read the data

channels = 1:2;

img = zeros([r.getSizeY() r.getSizeX() numel(channels) r.getSizeZ() r.getSizeT()], 'uint16');
for ti = 1:r.getSizeT()
    for cii = 1:numel(channels)
        for zi = 1:r.getSizeZ()
            img(:,:,cii,zi,ti) = bfGetPlane(r, r.getIndex(zi-1,channels(cii)-1,ti-1)+1);
        end
    end
end

%% read a FRAP region and visualize

ROIidx = 0; % omeMeta.getROICount() -> only one ROI for multiple FRAP
shapeIdx = 0; % but multiple shapes
x = {};
y = {};

for shapeIdx = 1:omeMeta.getShapeCount(ROIidx)/2
    
    shapeType = omeMeta.getShapeType(ROIidx,shapeIdx-1);
    if strcmp(shapeType, 'Polygon')
        
        disp(['reading polygon ' num2str(shapeIdx)]);
        
        s = strsplit(char(omeMeta.getPolygonPoints(ROIidx,shapeIdx-1)),{' ',','});
        s = cellfun(@str2double, s,'UniformOutput',false);
        s = cat(1, s{:});

        x{shapeIdx} = s(1:2:end);
        y{shapeIdx} = s(2:2:end);
        
    elseif strcmp(shapeType,'Ellipse')
        
        disp('ellipse?');
        X = omeMeta.getEllipseX(ROIidx, shapeIdx-1);
        Y = omeMeta.getEllipseY(ROIidx, shapeIdx-1);
        RX = omeMeta.getEllipseRadiusX(ROIidx, shapeIdx-1);
        RY = omeMeta.getEllipseRadiusY(ROIidx, shapeIdx-1);
        A = omeMeta.getEllipseTransform(ROIidx, shapeIdx-1);
        % this is a pain
    end
end


%% visualize

colors = lines(omeMeta.getShapeCount(ROIidx));

ti = 1;
zi = 1;
ci = S4Channel; %nucChannel;
imshow(img(:,:,ci,zi,ti),[1 2^12])
hold on
for shapeIdx = 1:omeMeta.getShapeCount(ROIidx)/2
    plot(x{shapeIdx},y{shapeIdx},'Color',colors(shapeIdx,:),'LineWidth',2)
end
hold off

saveas(gcf, fullfile(dataDir, [barefname '_FRAPframe' num2str(ti) '.png']));

%% read out profile in the mask

frapframe = 4; % actually 4

valnorm = {};

figure,
hold on 

for shapeIdx = 1:omeMeta.getShapeCount(ROIidx)/2
    
    mask = poly2mask(x{shapeIdx}, y{shapeIdx}, r.getSizeX(), r.getSizeY());
    mask = imerode(mask, strel('disk', 30));
    
    ci = S4Channel;

    val = zeros([1 r.getSizeT()]);
    for ti = 1:r.getSizeT()
        im = img(:,:,ci,zi,ti);
        val(ti) = mean(im(mask));
    end

    t = (1:r.getSizeT())*tres;
    valnorm{shapeIdx} = (val - min(val))/(max(val) - min(val));

    plot(t(frapframe:end),valnorm{shapeIdx}(frapframe:end))
end

hold off
xlim([0 t(end)]);
xlabel('time (sec)');
ylabel('recovered fraction')

%%
ti = 400;%r.getSizeT(); zi = 1; 
I = mat2gray(img(:,:,ci,zi,ti));
figure, imshow(cat(3,I, I + mask - imerode(mask,strel('disk',3)), I));

%% fitting

% function f to be fitted
%
% f = R*(1 - exp(-k t))
%
% R recovery fraction
% p = [R k]

tmax = {360,360};%{r.getSizeT(), r.getSizeT()};

p = {};

for shapeIdx = 1:omeMeta.getShapeCount(ROIidx)/2
    
    f = @(p,t) p(1)*(1 - exp(-p(2)*t));
    E = @(p) f(p,t(frapframe:tmax{shapeIdx})) - valnorm{shapeIdx}(frapframe:tmax{shapeIdx});

    R0 = 0.5;
    k0 = 0.01;

    pinit = [R0 k0];

    % lower and upper bounds
    lb(1) = 0;
    ub(1) = 1;
    lb(2) = 0;
    ub(2) = Inf;

    options.Algorithm = 'trust-region-reflective';
    options.TolFun = 1e-1;
    options.Display = 'off';

    [p{shapeIdx},resnorm,~,~,~,~,J] = lsqnonlin(E, pinit, lb, ub, options);
    sum(E(p{shapeIdx}).^2)
%     [Q,R] = qr(J,0);
%     Rinv = inv(R);
%     se = sqrt(sum(Rinv.*Rinv,2)*resnorm/dfe);

    pcur = p{shapeIdx};
    disp(['recovery fraction: ' num2str(pcur(1),2)]);
    tau=1/(60*pcur(2));
    disp(['kin + kout: ' num2str(pcur(2),2) ' s^-1 --> tau = ' num2str(tau,2) ' min']);
end

%%
clf

lw = 2;
hold on 

for shapeIdx = 1:omeMeta.getShapeCount(ROIidx)/2
    plot(t(frapframe:end),valnorm{shapeIdx}(frapframe:end),'Color',colors(shapeIdx,:),'LineWidth',lw)
end

legendstr = {};
%pplot = [0.25 0.035];
for shapeIdx = 1:omeMeta.getShapeCount(ROIidx)/2
    pplot = p{shapeIdx};
    plot(t(frapframe:tmax{shapeIdx}),f(pplot,t(frapframe:tmax{shapeIdx})),'Color',colors(shapeIdx,:),'LineWidth',lw)
    tau = 1/(60*p{shapeIdx}(2));
    legendstr{shapeIdx} = ['A = ' num2str(p{shapeIdx}(1),2) ', \tau=' num2str(tau,2) ' min'];
end

hold off
xlim([0 t(end)]);
ylim([0 0.7]);
legend(legendstr, 'Location','SouthEast');

fs = 15;
xlabel('time (sec)', 'FontSize',fs, 'FontWeight','Bold');
ylabel('recovered fraction N_{FRAP}/N', 'FontSize',fs, 'FontWeight','Bold');
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
title(barefname); 

saveas(gcf,fullfile(dataDir, [barefname '_FRAPplot']));
saveas(gcf,fullfile(dataDir, [barefname '_FRAPplot.png']));

%%

r.close();
