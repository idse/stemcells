function [a,di] = iwishrnd(sigma,df,di)
%IWISHRND Generate inverse Wishart random matrix
%   W=IWISHRND(SIGMA,DF) generates a random matrix W from the inverse
%   Wishart distribution with parameters SIGMA and DF.  The inverse of W
%   has the Wishart distribution with covariance matrix inv(SIGMA) and DF
%   degrees of freedom.
%
%   W=IWISHRND(SIGMA,DF,DI) expects DI to be lower triangular so that
%   DI'*DI = INV(SIGMA), i.e., the transpose of inverse of the Cholesky
%   factor of SIGMA. If you call IWISHRND multiple times using the same
%   value of SIGMA, it's more efficient to supply DI instead of computing
%   it each time.
%
%   [W,DI]=IWISHRND(SIGMA,DF) returns DI so it can be used again in
%   future calls to IWISHRND.
%
%   Note that different sources use different parameterizations for the
%   inverse Wishart distribution.  This function defines the parameter
%   SIGMA so that the mean of the output matrix is SIGMA/(DF-K-1), where
%   K is the number of rows and columns in SIGMA.
%
%   See also WISHRND.

%   Copyright 1993-2006 The MathWorks, Inc. 
%   $Revision: 1.4.4.9 $  $Date: 2007/05/23 19:15:37 $

% Error checking
if nargin<2
   error('stats:iwishrnd:TooFewInputs','Two arguments are required.');
end

[n,m] = size(sigma);
if n~=m
   error('stats:iwishrnd:BadCovariance','Covariance matrix must be square.');
end
if (~isscalar(df)) || (df<=0)
   error('stats:iwishrnd:BadDf',...
         'Degrees of freedom must be a positive scalar.')
end
if (df<n) % require this to ensure invertibility
   error('stats:iwishrnd:BadDf',...
         'Degrees of freedom must be no smaller than the dimension of SIGMA.');
end

% Get Cholesky factor for inv(sigma) unless that's already done
if nargin<3
   [d,p] = cholcov(sigma,0);
   if p~=0
      error('stats:iwishrnd:BadCovariance',...
            'Covariance matrix must be symmetric and positive definite.');
   end
   di = d'\eye(size(d));
else
   if ~isequal(size(di),size(sigma))
      error('stats:iwishrnd:BadCovFactor',...
            'DI must be the same size as SIGMA.')
   end
end

% Note:  the following would be more correct using inv(sigma) as the
% first argument, but the first argument will not be used if we pass
% in the factor di as the third argument.
a = wishrnd(sigma,df,di);
a = a\eye(size(a));
