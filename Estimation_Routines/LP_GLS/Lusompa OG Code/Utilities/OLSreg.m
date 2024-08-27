
function OLSresults = OLSreg(y,x,determ)

% Regresses input y on input x using OLS.
%
% Usage:   OLSresults = OLSreg(y,x,determ)
%
% Inputs:  y:      Vector of dependent variable    (T x 1)
%          x:      Matrix of independent variables (T x N)
%          determ: (opt.) Deterministic regressors: empty=none, 1=const,
%                  2=const&trend, 3=const&trend^2 (default is 0)
%
% Outputs: Structure "OLSresults" containing:
%   ~.y:     dependent vector (T x 1)
%   ~.x:     regressor matrix (T x N)
%   ~.nobs:  # observations
%   ~.nvar:  # variables
%   ~.beta:  bhat (N x 1)
%   ~.yhat:  yhat (T x 1)
%   ~.resid: residuals (T x 1)
%   ~.SSR:   sum of squared residuals
%   ~.bstd:  std deviations for bhat (N x 1)
%   ~.bint:  vector with 95% conf. intervals on beta (N x 2)
%   ~.tstat: t-stats  (N x 1)
%   ~.Rsqr:  R^2 (scalar)
%   ~.Rbar:  R_bar^2 (scalar)
%   ~.DW:    Durbin-Watson statistic
%   ~.F:     F-test for inclusion of determ into regression
%
% by B. Kolb, Jan. 2015, based on code by J.P. LeSage and B. Dillon (ols.m) 
%                        and A. Cesa Bianchi (OLSmodel.m)
% requires A. Holtsberg's tdis_inv.m, beta_inv.m and beta_pdf.m


%% Default determ
if ~exist('determ','var')
    determ = 0;
end

%% Add deterministic terms to regressors
switch determ
    case 0 
    case 1 % constant
        const = ones(size(y,1),1);
        x = [const x];
    case 2 % time trend and constant
        const = ones(size(y,1),1);
        trend = (1:size(x,1))';
        x = [const trend x];
    case 3 % squared time trend, linear time trend, and constant
        const = ones(size(y,1),1);
        trend = (1:size(x,1))';
        x = [const trend trend.^2 x];
    otherwise % not valid
        disp('Please select between 0 and 3 deterministic terms')
end

%% Get dimensions
if isempty(x)
    error('No regressors specified')
else
    [T, N] = size(x); 
    if (T ~= size(y,1)); error('x and y must have same # obs'); end
end

%% Inverse of (x'*x)
if T < 10000
  [~, r] = qr(x,0);
  xx_inv = (r'*r)\eye(N);
else
  xx_inv = (x'*x)\eye(N); % note: (x'*x)\eye(N) faster, but not exact!
end;

%% Get statistics and store them
bet   = xx_inv*(x'*y);
resid = y - x*bet;
SSR   = (resid'*resid)/(T-N);
SST   = (y - mean(y))'*(y - mean(y))/(T-1);
stdvb = sqrt(SSR*diag(xx_inv));
tcrit = -tdis_inv(.025,T);
resdf = resid(2:T) - resid(1:T-1);

OLSresults.y     = y;
OLSresults.x     = x;
OLSresults.nobs  = T;
OLSresults.nvar  = N;
OLSresults.beta  = bet;
OLSresults.yhat  = x*bet;
OLSresults.resid = resid;
OLSresults.SSR   = SSR;
OLSresults.bstd  = stdvb;
OLSresults.bint  = [OLSresults.beta-tcrit.*stdvb, ...
                    OLSresults.beta+tcrit.*stdvb];
OLSresults.tstat = OLSresults.beta./stdvb;
OLSresults.Rsqr  = 1 - (T-N)/(T-1)*SSR/SST;
OLSresults.Rbar  = 1 - SSR/SST;
OLSresults.DW    = (resdf'*resdf)/SSR/(T-N);

%% F-test
if determ > 0
    fx      = x(:,1); 
    fxx_inv = (fx'*fx)\1;
    fbeta   = fxx_inv*(fx'*y);
    fyhat   = fx*fbeta;
    fresid  = y - fyhat;
    fSSR    = fresid'*fresid;
    fRsqr   = 1 - fSSR/SST;
    OLSresults.F = (T-N)/(1-N)*(fRsqr-OLSresults.Rsqr)/(1-OLSresults.Rsqr);
end

