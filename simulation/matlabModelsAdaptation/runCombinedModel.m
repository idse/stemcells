tmax = 10;
nstep = 10;
LtoUse = [0, 0.5, 1, 2, 5, 10];
abToUse = [0, 0.1, 1, 5, 10];
n = length(abToUse);

figure;
for ii = 1:n
    for jj = 1:n
        subplot(n,n,(ii-1)*n+jj);
        doseResponse(tmax,nstep,abToUse(ii),abToUse(jj),LtoUse);
        title(['a = ' num2str(abToUse(ii)) '; b = ' num2str(abToUse(jj))]);
        xlim([0 15]);
    end
end
export_fig('doseresponse.eps','-transparent');

ligandlevel = 10;
figure;
for ii = 1:n
    for jj = 1:n
        subplot(n,n,(ii-1)*n+jj);
        stepVsRamp(ligandlevel,tmax,nstep,abToUse(ii),abToUse(jj));
        title(['a = ' num2str(abToUse(ii)) '; b = ' num2str(abToUse(jj))]);
    end
end
export_fig('stepvramp.eps','-transparent');


function stepVsRamp(ligandLevel,tmax,nstep,a,b)
modToRun = @(t,x) simpleCombinedStep(t,x,ligandLevel,tmax,nstep,a,b);
solStep = ode23(modToRun,[0 tmax*1.5],[0 0]);
modToRun = @(t,x) simpleCombinedRamp(t,x,ligandLevel,tmax,nstep,a,b);
solRamp = ode23(modToRun,[0 tmax*1.5],[0 0]);

 hold on;
plot(solStep.x,solStep.y(1,:),'LineWidth',3);
plot(solRamp.x,solRamp.y(1,:),'LineWidth',3);
xlim([0 1.5*tmax]);
end


function doseResponse(tmax,nstep,a,b,LtoUse)
hold on;
nL = length(LtoUse);
for ii = 1:nL
    modToRun = @(t,x) simpleCombinedStep(t,x,LtoUse(ii),tmax,nstep,a,b);
    sol = ode23(modToRun,[0 tmax*1.5],[0 0]);
    plot(sol.x,sol.y(1,:),'LineWidth',3);
    leg{ii} = ['L = ' num2str(LtoUse(ii))];
end
%legend(leg);
hold off;
end


function dx = simpleCombinedStep(t,x,Lfinal,tmax,nstep,a,b)
alpha = 1;
delta = 20;
beta = 1;
K = 2;
d = 0.1;

if t < tmax/2
    L = 0;
else
    L = Lfinal;
end

dx = zeros(2,1);
dx(1) = alpha - delta*x(1)/(K+L)-beta*x(1)*x(2);
dx(2) = a*L/(K+L)+b*x(1)-d*x(2);
end

function dx = simpleCombinedRamp(t,x,Lfinal,tmax,nstep,a,b)
alpha = 1;
delta = 20;
beta = 1;
K = 2;
d = 0.1;

if t < tmax/2
    L = 0;
elseif t < tmax
    t2 = t - tmax/2+tmax/(2*nstep);
    L = Lfinal*(t2-mod(t2,tmax/(2*nstep)))/(tmax/2);
else
    L = Lfinal;
end

dx = zeros(2,1);
dx(1) = alpha - delta*x(1)/(K+L)-beta*x(1)*x(2);
dx(2) = a*L/(K+L)+b*x(1)-d*x(2);
end