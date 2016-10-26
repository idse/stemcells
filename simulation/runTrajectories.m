function [X, t, Rfrac] = runTrajectories(params,b)

Ncells = 500;
dt = 0.01;
tmax = 5;
timax = ceil(tmax/dt);
D = 0.05; % noise strength;
c = 1;
t = dt*(1:timax);

npoly = length(params);
fparams = zeros(length(params)-1,1);
for ii = 1:length(params)-1
    fparams(ii) = (npoly-ii)*params(ii);
end

F = @(x,b) -polyval(fparams,x)-b;

x0 = zeros([Ncells 1]);
noise = randn([Ncells timax]);
X = zeros([Ncells timax]);

X(:,1) = x0;
for ti = 2:timax
    X(:,ti) = X(:,ti-1) + c*F(X(:,ti-1),b)*dt + sqrt(2*D*dt)*noise(:,ti); 
end

Rfrac = sum(X(:,end) < 0)/Ncells;
