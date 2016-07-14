% Plots a line graph showing pulse concentration.
% 
% logDir is the directory of the log file
% logName is the name of the log file
% startTime is the time of day in minutes that the experiment began. This
%       controls the offset of the pulses.
% low is the value of blank media
% high is the value of a pulse
function plotPulse(logDir,logName,startTime,low,high)

% read the time data
time = textread(fullfile(logDir,logName),'%d','delimiter',',');
time = (time - time(1))/60; % normalize to zero and convert to minutes

% populate the concentration vector
pulse = ~mod(round((1:numel(time))/2),2); % time points of pulses
level = zeros(numel(time),1);             % set up vector
level(pulse) = high;                      % pulses get high value
level(~pulse) = low;                      % blanks get low value

% plotting
plot(time-startTime,level,'LineWidth',2);

return