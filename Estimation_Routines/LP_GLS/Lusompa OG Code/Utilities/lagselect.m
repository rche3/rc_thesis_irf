
function [qAIC, qHQC, qBIC, qLBQ, qLLR]= lagselect(y,maxlag,determ,outp)

% Calculates lag selection criteria for VARs.
% [qAIC, qHQC, qBIC, qLBQ, qLLR]= lagselect(y,maxlag,outp)
% Inputs:  y:      variables
%          maxlag: maximum lag to consider
%          determ: (opt.) Deterministic regressors: 0=none, 1=const,
%                  2=const&trend, 3=const&trend^2 (default is 1)
%          outp:   (opt.) Result table is plotted if outp = 1
% Outputs: AIC: lag length proposed by Akaike Information Criterion (IC)
%                (Akaike 1973)
%          HQC: lag length proposed by Hannan/Quinn IC (Hannan/Quinn 1979)
%          BIC: lag length proposed by Bayesian IC (Schwarz 1979)
%          LBQ: lowest lag length for which zero cross-correlation of
%               residuals is no longer rejected in a Ljung-Box Q test
%          LLR: lowest lag length for which the null of a likelihood ratio
%               test (no difference to unrestricted) is rejected
% For details see Canova (2007), pp. 118 ff.
% by B. Kolb, Feb. 2015, based partially on code by F. Canova

%% Housekeeping
if size(y,1) < size(y,2)
    y = y';
end

% set defaults
if     nargin == 2
    determ = 1;
    outp = 1;
elseif nargin == 3
    outp = 1;
elseif nargin >= 5
    disp('Cannot handle more than 4 inputs')
end

%% Preallocation
AIC = zeros(maxlag,1); % Akaike Info. Crit.
HQC = zeros(maxlag,1); % Hannan & Quinn Info. Crit.
BIC = zeros(maxlag,1); % Bayesian Info. Crit.
LBQ = zeros(maxlag,1); % Ljung-Box Q test p-value
LLR = zeros(maxlag,1); % Likelihood ratio test p-value
nrg = zeros(maxlag,1);
df  = zeros(maxlag,1);
L   = zeros(maxlag,1);

LLR(maxlag) = NaN; % "end-point problem"

%% Calculate statistics
for ii = 1:maxlag
    
    VARresults = VARreg(y,ii,determ);

    T     = VARresults.nobs;
    N     = VARresults.nvar;
    sigma = VARresults.sigma;
    u     = VARresults.resid;
    A     = VARresults.Ft;
    
    % AIC, HQC, BIC
    AIC(ii) = log(det(sigma)) + 2*N*(ii*N+determ)/T;
    HQC(ii) = log(det(sigma)) + 2*log(log(T))*N*(ii*N+determ)/T;
    BIC(ii) = log(det(sigma)) + log(T)*N*(ii*N+determ)/T;
    
    % Ljung-Box Q statistic
    Q1    = 0;
    for jj = 1:maxlag
        uj = [zeros(jj,size(u,2)); u(1:end-jj,:)];
        Cj = (u'*uj)/T;
        Q1 = Q1 + trace(Cj'/sigma*Cj/sigma) / (T-jj);
    end
    
    Q1      = (T^2)*Q1;
    dfQ     = N^2*(maxlag-ii);
    LBQ(ii) = 1 - chi2cdf(Q1,dfQ);
    
    % Likelihood ratio test
    L(ii)       = -(T*N/2)*(1+log(2*pi)) + (T/2)*log(det(sigma^(-1)));
    if ii >= 2
    LLR(ii-1) = cdf('chi2',2*(T-N*ii)/T*(L(ii)-L(ii-1)),N);
    end
    
    % for displaying
    nrg(ii) = size(A,1) - 1;
    df(ii)  = T - size(A,1); 
    
end

%% Selected number of lags by each criterion
[~,qAIC] = min(AIC);
[~,qHQC] = min(HQC);
[~,qBIC] = min(BIC);
qLBQ     = find(LBQ<=0.05,1);
qLLR     = find(LLR>=0.05,1) - 1;


%% Display Output
if outp == 1
    q = (1:1:maxlag)';
    disp('                                       Selecting the Lag Order ');
    disp(' -------------------------------------------------------------------------------------------------');
    disp('     Lag       AIC       HQC       BIC    LBQ(pval) LLR(pval)  log(L)      nrg       df');
    disp([q AIC HQC BIC LBQ LLR L nrg df]);
    disp(' ');
    disp('    Selected Number of Lags');
    disp(' ----------------------------------');
    disp('    AIC   HQC   BIC   LBQ   LLR');
    disp([qAIC qHQC qBIC qLBQ qLLR]);
end

