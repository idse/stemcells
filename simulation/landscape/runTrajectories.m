function [X, t, Rfrac] = runTrajectories(params,b,D)
% Function to run trajectories in a 1D landscape defined by a polynomial
% Inputs:
% -params defines a polynomial as for polyval of order length(params)-1 with
% the first entry being the highest order coefficient. 
% -b is the magnitude of the constant force
% -D is the magnitude of the stochastic noise
% Outputs: 
% X - trajectories - each row is a cell and each column a different time
% t - the times corresponding to the columns of X
% Rfrac - the fraction that end up in the left basin. 
%
% see also polyval, plotTrajectories

if ~exist('b','var');
    b = 0;
end

if ~exist('D','var')
    D = 0.05;
end

Ncells = 100;
dt = 0.01;
tmax = 5;
timax = ceil(tmax/dt);
c = 1;
t = dt*(1:timax);

%Take the derivative of the landscape
npoly = length(params);
fparams = zeros(length(params)-1,1);
for ii = 1:length(params)-1
    fparams(ii) = (npoly-ii)*params(ii);
end

%define the force function
F = @(x,b) -polyval(fparams,x)-b;

%do the intergration
x0 = zeros([Ncells 1]);
noise = randn([Ncells timax]);
X = zeros([Ncells timax]);

X(:,1) = x0;
for ti = 2:timax
    X(:,ti) = X(:,ti-1) + c*F(X(:,ti-1),b)*dt + sqrt(2*D*dt)*noise(:,ti); 
end

%Fraction of endo
Rfrac = sum(X(:,end) < 0)/Ncells;
