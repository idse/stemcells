function mmdKS_results = mmdKS(data1, data2, varargin)
% mmdKS.m
% Written February 27, 2011 by Marcello DiStasio. SUNY Downstate Medical
% Center.
% 
% mmdKS_results = mmdKS(data1, data2, [alpha])
%
%  Compute Kolmogorov-Smirnov test statistics for two data sets (1D row or column vectors).  
%  alpha is an optional confidence level (default alpha=0.05)
%  Computes confidence intervals around the CDF for data1 
%  (i.e. data1 is the true distribution to be matched by data2)
% 
%
%  The fields in the return struct are:
%    mmdKS_results.sampleCDF1 = empirical CDF of data1
%    mmdKS_results.sampleCDF2 = empirical CDF of data1
%    mmdKS_results.D          = K-S statistic (max difference between CDFs);
%    mmdKS_results.H          = Reject or accept the null hypothesis (0 = same distribution)   
%    mmdKS_results.pValue     = pValue (probability of these two CDFs being the same)
%    mmdKS_results.CI_upper   = upper confidence bound (certainty = 1-alpha)
%    mmdKS_results.CI_lower   = lower confidence bound (certainty = 1-alpha)

% Make sure the data are column vectors:
if (size(data1,2) > 1)
    data1 = data1';
end
if (size(data2,2) > 1)
    data2 = data2';
end

% Table data for signficance.  Look for user supplied alpha value here
alphas = [0.10 0.05 0.025 0.01 0.005 0.001];
coeff_alphas = [1.22 1.36 1.48 1.63 1.73 1.95];

% Get the user specified alpha value's corresponding coefficient for
% critical value of the D statistic (c).  Default is alpha = 0.05.
if (nargin > 2)
    if any(alphas==varargin{1})
        sig_level = find(alphas==varargin{1});
    else
        fprintf('Supplied alpha value (%f) is not in the table.  Please choose one of:',varargin{3})
        disp(alphas)
        return
    end
else
    sig_level = 2;
end
alpha = alphas(sig_level);
c     = coeff_alphas(sig_level);

tail  =  0; % Two-sided test

% Calculate F1(x) and F2(x), the empirical (i.e., sample) CDFs.
binEdges    = [-inf ; sort([data1;data2]) ; inf];

binCounts1  = histc (data1, binEdges, 1);
binCounts2  = histc (data2, binEdges, 1);

sumCounts1  = cumsum(binCounts1)./sum(binCounts1);
sumCounts2  = cumsum(binCounts2)./sum(binCounts2);

sampleCDF1  = sumCounts1(1:end-1);
sampleCDF2  = sumCounts2(1:end-1);

% Compute the test statistic of interest.
switch tail
   case  0      %  2-sided test: T = max|F1(x) - F2(x)|.
      deltaCDF  =  abs(sampleCDF1 - sampleCDF2);

   case -1      %  1-sided test: T = max[F2(x) - F1(x)].
      deltaCDF  =  sampleCDF2 - sampleCDF1;

   case  1      %  1-sided test: T = max[F1(x) - F2(x)].
      deltaCDF  =  sampleCDF1 - sampleCDF2;
end

KSstatistic   =  max(deltaCDF);


% Compute confidence measures:
n1     =  length(data1);
n2     =  length(data2);
n      =  n1 * n2 /(n1 + n2);

D_critical = c*sqrt(1/n);

% Compute the asymptotic P-value approximation and accept or
% reject the null hypothesis on the basis of the P-value.
%
lambda =  max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * KSstatistic , 0);

if tail ~= 0        % 1-sided test.
   pValue  =  exp(-2 * lambda * lambda);
else
    % 2-sided test (default).
    %
    %  Use the asymptotic Q-function to approximate the 2-sided P-value.
    %
   j       =  (1:101)';
   pValue  =  2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
   pValue  =  min(max(pValue, 0), 1);

end

% If alpha exceeds the pvalue, then report the result as significant
% i.e. H = 1 means REJECT THE NULL HYPOTHESIS 
% H_0 = Datasets are drawn from DIFFERENT POPULATIONS
% That is, H = 1 means that these two datasets are drawn from the same
% distribution.
H  =  (alpha >= pValue);



% Generate return struct
mmdKS_results.sampleCDF1 = sampleCDF1;
mmdKS_results.sampleCDF2 = sampleCDF2;
mmdKS_results.D          = KSstatistic;
mmdKS_results.H          = H;   
mmdKS_results.pValue     = pValue;
mmdKS_results.CI_upper   = sampleCDF1 + D_critical;
mmdKS_results.CI_lower   = sampleCDF1 - D_critical;

% 
% % Plot results
% figure();
% plot(mmdKS_results.sampleCDF1); hold on; 
% plot(mmdKS_results.sampleCDF2,'k'); 
% plot(mmdKS_results.CI_upper,'r'); 
% plot(mmdKS_results.CI_lower,'r');
% title(sprintf('Kolmogorov-Smirnov test. D = %f, p-value = %f',mmdKS_results.D, mmdKS_results.pValue))

end