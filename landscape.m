% first attempt at landscape simulation, 1D

% potential:
% V = x^4/4 + a x^2/2   
% a is the landscape parameter

x = -2:0.1:2;
a = 2;
b = 0;
V = @(x) x.^4/4 - a*x.^2/2;

plot(x,V(x));

%% force
% F = -dV/dx + b

F = @(x,b) -x.^3 + a*x - b;

plot(x,V(x) + b.*x,'r')
hold on
plot(x,F(x,b))
hold off
legend('potential','force')

%% Euler-Maruyama integration
%
% dx/dt = c F
% we assume we can place the initial condition in the middle of the
% symmetric background landscape

b = 0;
c = 1;
Ncells = 100;
dt = 0.01;
tmax = 5;
timax = ceil(tmax/dt);
t = dt*(1:timax);
D = 0.05; % noise strength;

x0 = zeros([Ncells 1]);
noise = randn([Ncells timax]);
X = zeros([Ncells timax]);

X(:,1) = x0;
for ti = 2:timax
    X(:,ti) = X(:,ti-1) + c*F(X(:,ti-1),b)*dt + sqrt(2*D*dt)*noise(:,ti); 
end

% visualize trajectories
plot(x,V(x)+b.*x)
hold on
plot(X',t)
hold off
xlim([x(1) x(end)])

%% visualize final distribution

hist(X(:,end),50)

sum(X(:,end)>0)
sum(X(:,end)<0)

%% fraction in right well vs b

N = 15;
bval = linspace(-1,1,N);
Rfrac = zeros([1 N]);

for bi = 1:N;

    b = bval(bi);
    
    x0 = zeros([Ncells 1]);
    noise = randn([Ncells timax]);
    X = zeros([Ncells timax]);

    X(:,1) = x0;
    for ti = 2:timax
        X(:,ti) = X(:,ti-1) + c*F(X(:,ti-1),b)*dt + sqrt(2*D*dt)*noise(:,ti); 
    end
    
    Rfrac(bi) = sum(X(:,end) > 0)/Ncells;
end

plot(bval,Rfrac)
