% runs original model of IH with new functions
params=  [0.25 0 -1 0 0];
b = 0;
[X, t, Rfrac]=runTrajectories(params,b);
plotTrajectories(X,t,params,b);
%% produce force param plot
N = 25;
bval = linspace(-1,1,N);
Rfrac = zeros([1 N]);

for bi = 1:N;
    b = bval(bi);
    [~,~,Rfrac(bi)]=runTrajectories(params,b);
end
plot(bval,Rfrac)
%% compare data and simulation
act = [0 0.5 1 5 10 50 100];
endo = [0 0 3 32 38 82 88];
N = length(endo);
Rfrac = zeros([1 N]);
bval = -1+(act/100*2);
N = length(bval);
for bi = 1:N;
    b = bval(bi);
    [~,~,Rfrac(bi)]=runTrajectories(params,b);
end
plot(bval,Rfrac,'.-','LineWidth',3,'MarkerSize',16); hold on;
plot(bval,endo/100,'.-','LineWidth',3,'MarkerSize',16); hold off;
%% multiparam fit
aval = 1;
params = [0.25 0 -aval 0 0];
for bi = 1:N;
    b = bval(bi);
    [~,~,Rfrac(bi)]=runTrajectories(params,b);
end
datadiff = Rfrac-endo/100;
datadiff=sum(datadiff.^2);
Ntrial = 1000;
Temp = 0.01; %little bit of thermal noise can get fitting unstuck
steps_since_update = 0;
bestval = datadiff;
bestparams = params;
D = 0.05;
bestD = D;
npoly = length(params);

paramsToFit = [1 3 4 6];
for ii = 1:Ntrial
    
    %pick the parameter to change and adjust it
    r0 = randi(length(paramsToFit));
    paramToChange = paramsToFit(r0);
    r1 = 0.1*(2*rand()-1); %random between -0.1 and 0.1
    
    if paramToChange > npoly
        Dnew = D + r1*0.1; %change noise magnitude
        paramsnew = params;
    else
        paramsnew = params;
        Dnew = D;
        paramsnew(paramToChange) = params(paramToChange)+r1;
    end
    %get the predicted values
    for bi = 1:N;
        b = bval(bi);
        [~,~,Rfrac(bi)]=runTrajectories(paramsnew,b,Dnew);
    end
    
    %compare to data and decide whether to accept
    datadiff2 = Rfrac-endo/100;
    datadiff2=sum(datadiff2.^2);
    if datadiff2 < datadiff || 0.5*exp(-(datadiff2-datadiff)/Temp) > rand()
        D = Dnew;
        params = paramsnew;
        datadiff = datadiff2;
    end
    if datadiff < bestval
        bestval = datadiff;
        bestparams = params;
        bestD = D;
        steps_since_update = 0;
    else
        step_since_update = steps_since_update + 1;
    end
    
    if steps_since_update > 200
        break;
    end
    disp(['Trial ' int2str(ii) ' aval = ' num2str(aval), ' datadiff = ' num2str(datadiff)]);
end
%% compare data and simulation
act = [0 0.5 1 5 10 50 100];
endo = [0 0 3 32 38 82 88];

N = length(endo);
Rfrac = zeros([1 N]);
bval = -1+(act/100*2);
N = length(bval);
for bi = 1:N;
    b = bval(bi);
    [~,~,Rfrac(bi)]=runTrajectories(bestparams,b,bestD);
end
plot(bval,Rfrac,'.-','LineWidth',3,'MarkerSize',16); hold on;
plot(bval,endo/100,'.-','LineWidth',3,'MarkerSize',16); hold off%% show simulation of the final landscape, no bias
%% best landscape, 0 force
[X, t]=runTrajectories(bestparams,0,bestD);
plotTrajectories(X,t,bestparams,0);
%%
figure;
plot(bval,Rfrac,'.-','LineWidth',3,'MarkerSize',16); hold on;
plot(bval,endo/100,'.-','LineWidth',3,'MarkerSize',16); hold off;
legend({'Model','Data'},'FontSize',24);


