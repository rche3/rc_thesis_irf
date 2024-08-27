function [betahat,vbeta,se_beta,ser,rbarsq,q_df] = har_ols(y,x,q);
%{
     Modified by MWW, 6-10-2018

     Procedure for estimating the regression y = xbeta+ u
     The procedure produces the OLS estimate of b
     and a hetero/autocorrelation robust estimate of
     the covariance matrix for betahat using Cosine transforms as basis function
     (Mueller (2004)). 

Input:
     y = Tx1
     x = Txk
     q = number of df (number of cosine terms used to compute covariance matrix)
         (when q = 0, use the default value q = 0.41*T^(2/3) 

Output:
    Beta = OLS estimate of beta (kx1)
    VBeta = Robust estimate of covariance matrix of beta (kxk)
    se_beta = standard errors of betahat (note: (betahat-beta_true)/se_beta ~ Student-t with q d.f.
    ser = standard error of regression
    rbarsq = adjusted R-squared
    q_df = degrees of freedom for HAR estimate of covariance matrix. Note: q_dr = input value q unless q is computed internally

%}
xx=x'*x;
xxi = inv(xx); 
betahat=xxi*x'*y;
uhat=y-x*betahat;
z = x.*repmat(uhat,1,size(x,2));    % X*uhat
T = size(y,1);

if q == 0;
  q = floor(0.41*T^(2/3));
end

% Cosine weights
fvec = pi*(1:1:q);
tvec = (0.5:1:T)'/T;
psi = (sqrt(2)/sqrt(T))*cos(tvec*fvec);  % normalized so that psi'psi = 1

% Construct Estimate of Omega
w = psi'*z;
omega_hat_mat = w'*w/q;

vbeta=T*xxi*omega_hat_mat*xxi';

% Compute Other Statistics of interest
y_m = y - mean(y);
tss = y_m'*y_m;
ess = uhat'*uhat;
ndf = size(x,1)-size(x,2);
ser = sqrt(ess/ndf);
rbarsq = 1 - (ess/ndf)/(tss/(size(x,1)-1));
se_beta = sqrt(diag(vbeta));
q_df = q;

end