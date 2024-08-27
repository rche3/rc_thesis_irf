function [s0,s1,m] = nwbandwidth(xseries,eta_hat,start_IP,end_IP,nq)
% Note: This is a routine.
%
%
% See Newey and West (1994).
% This procedure is specialized for: OLS and single stochastic right hand side variable, which is the first
% lag of the x series.  Note the word "lag": x{1}(t)*etahat(t) is cross-produce of rhs and residual.
% In addtion to (1) xseries, variables passed are
% (2) eta_hat: series of least squares residuals
% (3) start_IP: start date of the regression
% (4) end_IP: end date of the regression
% (5) nq: known order of MA of the residual (and of cross-product of rhs variable and residual)
% 
% Variables return
% s0, s1: scalar real variables, as defined in Newey and West(1994).
% m: integer optimal bandwidth

Tnow=end_IP-start_IP+1; % number of data for the regression

x_lag=lagmatrix(xseries,1);
h_hat=x_lag(start_IP-nq:end_IP-nq).*eta_hat;

% compute mean square of h_hat (need to demean h_hat)
demean_h_hat=h_hat-mean(h_hat);
sigma0=(demean_h_hat'*demean_h_hat)/Tnow;
sigma=zeros(nq,1);

for j=1:nq
temp_h=h_hat(1+j:Tnow);
temp_hj_lagmatrix=lagmatrix(h_hat,j);
temp_hj=temp_hj_lagmatrix(1+j:Tnow);

demean_temp_h=temp_h-mean(temp_h);
demean_temp_hj=temp_hj-mean(temp_hj);
sigmanow=(demean_temp_h'*demean_temp_hj)/(Tnow-j) + mean(temp_h)*mean(temp_hj);
ratio=(Tnow-j)/Tnow;
sigma(j,1)=ratio*sigmanow;
clear temp_h temp_hj_lagmatrix temp_hj    
end

% compute s0 and s1
s0=sigma0+2*sum(sigma);
s1=0;
for j=1:nq
s1=s1+2*j*sigma(j);
end

% compute gamma
gamma=1.1147*((s1/s0)^2)^(1/3);

% compute m
real_m=gamma*(Tnow^(1/3));
m=floor(min(real_m,0.5*Tnow));


end









