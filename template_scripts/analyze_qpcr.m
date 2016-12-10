clear all; close all;

addpath(genpath('/Users/idse/repos/Warmflash/')); 
dataDir = '/Users/idse/data_tmp/qpcr/new_qpcr';
%dataDir = '/Users/idse/data_tmp/qpcr';

%filenames = {'161103 KB LEFTY NODAL SKI SKIL time series_data.txt'};
% filenames = {   '161101 KB SKIL Time Series_3_data.txt',...
%                 '161103 KB LEFTY NODAL SKI SKIL time series_data.txt'};

filenames = {'161202 KB ATP5O PBGD LEFTY1 NODAL1 time series comparison_data.txt'};

filenames = {'161205 KBPBGD ATP5O SNAI1 PAI1 Activin time series comparison_data.txt',...
            '161205 KB SMAD7 MIXL1 SKI SKILActiving time series comparison_data.txt'};

%%
% SOP : StepOne Plus
data = {};
for i = 1:numel(filenames)
    data{i} = readSOPdata3(fullfile(dataDir,filenames{i}));
end

%% collect CT values

CTmean = [];
CTstd = [];
Ntargets = 0;
Nsamples = 0;
targets = [];

for i = 1:numel(data)
	
    D = data{i};
    
	CTmean = cat(2, CTmean, D.CT);
    CTstd = cat(2, CTstd, D.CTstd);
    
    Ntargets = Ntargets + D.Ntargets;
    targets = cat(1, targets, D.targets);
    
    if i == 1
        Nsamples = D.Nsamples;
    elseif Nsamples ~= D.Nsamples
        error('samples dont match');
    end
end

%%
refIdx = find(strcmp(targets,'ATP5O'));
CTref = CTmean(:,refIdx);

CTnorm = zeros([Nsamples Ntargets]);
for targeti = 1:Ntargets
	CTnorm(:,targeti) = -(CTmean(:,targeti) - CTref);
    CTnorm(:,targeti) = CTnorm(:,targeti) - CTnorm(1,targeti);
end

%% raw value

times = [0 1 2 4 8 12 24];
%times = [0 2 4 8];

clf 
bar(times, CTmean)
legend(targets,'Location','SouthEast')
ylim([20 35]);

%% normalized

% bar(times, CTnorm(:,2:end))
% legend(D.targets(2:end))

targeti = setdiff(1:Ntargets, refIdx);

clf
color = lines(numel(targeti));
hold on
for i = 1:numel(targeti)
    plot(times, CTnorm(:,targeti(i)),...
       '-x','Color',color(i,:),'LineWidth',2);
    %errorbar(times, CTnorm(:,targeti(i)), CTstd(:,targeti(i)),...
    %    '-x','Color',color(i,:),'LineWidth',2);
end
hold off
xlim([0 24]);
ylim([-0.5 3.5]);
legend(targets(targeti),'location','NorthEast')

fs = 15;
xlabel('time (hrs)', 'FontSize',fs, 'FontWeight','Bold');
ylabel('C_T - C_{T,ATP5O}', 'FontSize',fs, 'FontWeight','Bold');
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')

%saveas(gcf, fullfile(dataDir, 'relCT');

%% normalized exponential

clf
targeti = setdiff(1:Ntargets, refIdx);
color = lines(numel(targeti));
hold on
for i = 1:numel(targeti)
    X = CTnorm(:,targeti(i));
    plot(times, 2.^X,'-x','Color',color(i,:),'LineWidth',2)
end
hold off
legend(targets(targeti))


fs = 15;
xlabel('time (hrs)', 'FontSize',fs,'FontWeight','Bold');
ylabel('normalized level', 'FontSize',fs, 'FontWeight','Bold');
set(gcf,'color','w');
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
    