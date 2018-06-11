function [times, pulses] = pulseProfile(startFrame, endFrame, pulseLength, waitLength, dt)
    % create vectors for plotting pulse profile
    % 
    % [times, pulses] = pulseProfile(startFrame, pulseLength, waitLength, Cmax)
    %
    % use output: plot(times, Cmax*pulses)

    Npulses = numel(pulseLength);

    pulse = [0 0 1 1];
    pulses = [];
    times = [];
    t0 = -startFrame*dt;
    waitLength = [startFrame*dt waitLength(1:end-1)];
    for i = 1:Npulses

        t1 = t0 + waitLength(i);
        t2 = t1 + pulseLength(i);
        pulseTimes = [t0 t1 t1 t2];
        t0 = t2;
        times = [times pulseTimes];
        pulses = [pulses pulse];
    end
    times = [times times(end) endFrame*dt];
    pulses = [pulses 0 0];
    %pulses = (Cmax-Cmin)*pulses + Cmin;

    tmax = endFrame;
    tinit = (1 - startFrame)*dt;
    tfinal = (tmax - startFrame)*dt;
    t = ((1:tmax) - startFrame)*dt;

    %allpulset = [(0:Npulses-1)*tlength, (0:Npulses-1)*tlength + plength];
    %times = sort([tinit allpulset allpulset tfinal]);

end