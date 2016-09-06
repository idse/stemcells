% suppose we have a distribution of time traces
%
% one source of variation is translation and scaling relative to the mean
% close to the mean we can represent these as
% 
% (1 + a T)f and (1 + b D)f
% where f is the the time trace 
% and T and D are infinitesimal translation and dilation generators, i.e. 
% d/dt and t d/dt 
%
% if in principal component analysis we diagonalize the covariance
% matrix to decompose 
%
% f = sum_i c_i f_i 
%
% comparing to the above we see that f_i will be the mean, T acting on the
% mean (the time derivative of the mean), D acting on the mean, etc.
% 
% we can exclude these principal components and the remaining leading ones
% will give a basis for the more non-trivial variation of the time traces
% 
% it is then easy to check the distribution along those directions and
% cluster the trajectories based on that

%%
% generate a distribution of time traces randomly shifted relative to one
% another

nTime = 50;
nTraces = 100;
T = linspace(1,4,nTime);
timeTrace = zeros([nTraces nTime]);

for i = 1:nTraces
    a = 0.2*randn(1);
    timeTrace(i,:) = sin(2*pi*T+a);
end
    
plot(timeTrace')

%% do PCA

% nTime x nTime covariance matrix
S = cov(timeTrace);
[V,D] = eigs(S);

% the first two principal components
plot(V(:,1)')
hold on 
plot(V(:,2)')
hold off
axis([0 50 -1/4 1/4])

%% first 2 PC coincide perfectly with the mean and the derivative of the mean

% disagreement here is due to the one sided approximation of the derivative
% and edge effects of the convolution (and there is a sign difference)

% normalize to compare to eigenvector
meanTrace = mean(timeTrace)/sqrt(sum(mean(timeTrace).^2));
dMeanTrace = conv(meanTrace,[1 -1]/2);
dMeanTrace = -dMeanTrace/sqrt(sum(dMeanTrace.^2));

figure,
plot(meanTrace,'g')
hold on
plot(dMeanTrace,'g')
hold off
axis([0 50 -1/4 1/4])

