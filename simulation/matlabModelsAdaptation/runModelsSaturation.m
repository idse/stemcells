ligandLevel = 10;
tmax = 100;
nstep = 10;

modToRun = @(t,x) simpleFeedForwardStep(t,x,ligandLevel,tmax,nstep);
solFFStep = ode23(modToRun,[0 tmax*1.5],[0 0]);
modToRun = @(t,x) simpleFeedForwardRamp(t,x,ligandLevel,tmax,nstep);
solFFRamp = ode23(modToRun,[0 tmax*1.5],[0 0]);
modToRun = @(t,x) simpleFeedbackStep(t,x,ligandLevel,tmax,nstep);
solFBStep = ode23(modToRun,[0 tmax*1.5],[0 0]);
modToRun = @(t,x) simpleFeedbackRamp(t,x,ligandLevel,tmax,nstep);
solFBRamp = ode23(modToRun,[0 tmax*1.5],[0 0]);



figure; subplot(2,2,1); hold on;
plot(solFFStep.x,solFFStep.y(1,:),'LineWidth',3);
plot(solFFRamp.x,solFFRamp.y(1,:),'LineWidth',3);
xlim([0 1.5*tmax]);title('FF - y'); legend({'step','ramp'});

subplot(2,2,2); hold on;
plot(solFBStep.x,solFBStep.y(1,:),'LineWidth',3);
plot(solFBRamp.x,solFBRamp.y(1,:),'LineWidth',3);
xlim([0 1.5*tmax]); 
title('FB - y');

subplot(2,2,3); hold on;
plot(solFFStep.x,solFFStep.y(2,:),'LineWidth',3);
plot(solFFRamp.x,solFFRamp.y(2,:),'LineWidth',3);
xlim([0 1.5*tmax]); title('FF - x');

subplot(2,2,4); hold on;
plot(solFBStep.x,solFBStep.y(2,:),'LineWidth',3);
plot(solFBRamp.x,solFBRamp.y(2,:),'LineWidth',3);
xlim([0 1.5*tmax]); 
title('FB - x');

function dx = simpleFeedForwardStep(t,x,Lfinal,tmax,nstep)
a = 1; 
b = 2;
k = 0.1;
K = 2;
d = 0.5;
c = 20;

if t < tmax/2
    L = 0;
else
    L = Lfinal;
end

dx = zeros(2,1);
dx(1) = a - c*x(1)/(K+L)-b*x(1)*x(2);
dx(2) = d*L/(K+L) -k*x(2);
end

function dx = simpleFeedForwardRamp(t,x,Lfinal,tmax,nstep)
a = 1; 
b = 2;
k = 0.1;
K = 2;
d = 0.5;
c = 20;

if t < tmax/2
    L = 0;
elseif t < tmax
    t2 = t - tmax/2+tmax/(2*nstep);
    L = Lfinal*(t2-mod(t2,tmax/(2*nstep)))/(tmax/2);
else
    L = Lfinal;
end

dx = zeros(2,1);
dx(1) = a - c*x(1)/(K+L)-b*x(1)*x(2);
dx(2) = d*L/(K+L) -k*x(2);
end

function dx = simpleFeedbackStep(t,x,Lfinal,tmax,nstep)
a = 1; 
b = 2;
c = 20;
K = 2;
k = 0.1;
d = 10;

if t < tmax/2
    L = 0;
else
    L = Lfinal;
end

dx = zeros(2,1);
dx(1) = a - c*x(1)/(L+K)-b*x(1)*x(2);
dx(2) = d*x(1)-k*x(2);
end

function dx = simpleFeedbackRamp(t,x,Lfinal,tmax,nstep)
a = 1; 
b = 2;
c = 20;
K = 2;
k = 0.1;
d = 100;

if t < tmax/2
    L = 0;
elseif t < tmax
    t2 = t - tmax/2+tmax/(2*nstep);
    L = Lfinal*(t2-mod(t2,tmax/(2*nstep)))/(tmax/2);
else
    L = Lfinal;
end

dx = zeros(2,1);
dx(1) = a - c*x(1)/(L+K)-b*x(1)*x(2);
dx(2) = d*x(1)-k*x(2);

end