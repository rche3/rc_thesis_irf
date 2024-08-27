function [std_NW,t_NW] = NeweyWest1994(y,X,m)
% PURPOSE: computes OLS with Robust, Newey-West heteroscedastic-serial consistent standard errors
% modified OLS with Newey-West and Hansen-Hodrick SE from https://www.mathworks.com/matlabcentral/fileexchange/43259-ols-with-newey-west-and-hansen-hodrick-se


% Inputs:
%  y = T x 1 vector, left hand variable data
%  X = T x n matrix, right hand variable data
%  m = number of lags to include in NW corrected standard errors
%
%Note: you must make one column of X a vector of ones if you want aXX
%   constant.
% Output:
%  std_NW  = corrected standard errors.
%  t_NW   = t-stat for NW
%Note: For chi-square test program checks whether first is a constant and ignores that one for
%       test. If there is only one beta the program does not report X^2
%       results since t_stat^2= X2.
%Note: program automatically displays outputs in an organized format. If you want
%to disable the automatic display just comment lines 129-136.

%Estimate Betas and Residuals
[T,n]   =   size(X);
beta_YY    =   (inv(X'*X))*X'*y;
u       =   y-X*beta_YY;
u= u*ones(1,n);
err=X.*u; %estimating residuals for each beta



%Calculate NW Corrected Standard Errors
V=[err'*err]/T;
if m > -1
    for ind_i = (1:m);
        S = err(1:T-ind_i,:)'*err(1+ind_i:T,:)/T;
        V = V + (1-1*ind_i/(m+1))*(S + S');
    end;
end;

D       =   inv((X'*X)/T);
varb = 1/T*D*V*D;
seb = diag(varb);
std_NW = sign(seb).*(abs(seb).^0.5);

%NW
if X(:,1) == ones(size(X,1),1);
    chi2val = beta_YY(2:end,:)'*inv(varb(2:end,2:end))*beta_YY(2:end,:);
    dof = size(beta_YY(2:end,1),1);
    pval = 1-cdf('chi2',chi2val, dof);
    X2_NW(:,1:3) = [chi2val dof pval];
else;
    chi2val = beta_YY(:,:)'*inv(varb)*beta_YY(:,:);
    dof = size(beta_YY(:,1),1);
    pval = 1-cdf('chi2',chi2val, dof);
    X2_NW(:,1:3) = [chi2val dof pval];
end;

% T_Stats
t_NW =  beta_YY./std_NW;
clear V S D varb seb chi2val dof pval
%-----------------------------------------------------------------------%
%Calculate Robust Standard Errors
V=[err'*err]/T;
D       =   inv((X'*X)/T);
varb = 1/T*D*V*D;
seb = diag(varb);
std_R = sign(seb).*(abs(seb).^0.5);

%F test for all coeffs (except constant) zero -- actually chi2 test
%NW
if X(:,1) == ones(size(X,1),1);
    chi2val = beta_YY(2:end,:)'*inv(varb(2:end,2:end))*beta_YY(2:end,:);
    dof = size(beta_YY(2:end,1),1);
    pval = 1-cdf('chi2',chi2val, dof);
    X2_R(:,1:3) = [chi2val dof pval];
else;
    chi2val = beta_YY(:,:)'*inv(varb)*beta_YY(:,:);
    dof = size(beta_YY(:,1),1);
    pval = 1-cdf('chi2',chi2val, dof);
    X2_R(:,1:3) = [chi2val dof pval];
end;

% T_Stats
t_R =  beta_YY./std_R;

% Formatting
%disp('  beta_YY     se_NW     t_NW         ');
%disp([beta_YY std_NW  t_NW ]);
end
