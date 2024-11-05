clear all
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code')
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code/Monte Carlos')
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code/bias_code/procedures')
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code/Utilities')

%Generating data and true impulse responses for the bivariate VAR(1) used
%in the Additional Monte Carlo section of Online Appendix.

B1=[.9 0;.5 .5];

maxh = 15; %maximum horizon one wants to estimate for LP

q=size(B1,1);
trueirf=zeros(q,q,maxh);

for h=1:maxh
    trueirf(:,:,h)=(B1^h);
end


q=size(B1,1);


    T=250; %sample size
    burnin=1000; %burnin


Errors=randn(q,T+burnin);
y=zeros(q,T+burnin);
y(:,1)=randn(q,1);


for t=2:T+burnin
    y(:,t)=B1*y(:,t-1)+Errors(:,t);
end
y=y(:,burnin+1:end);

blocksize=1; %block size for the bootstrap
badj=1; %1 for bias adjustment 0 for no bias adjustment. The econometrics toolbox is required for bias adjustment. It's debatable whether you should bias adjustment see Efron and Tibshirani (1993). 
% So depending on what your doing you may or may not want to bias adjust. Uses bias adjustment from West and Zhao 2016
% see bias code folder and this link for more details https://www.ssc.wisc.edu/~kwest/publications/2010/Adjusting%20for%20Bias%20in%20Long%20Horizon%20Regressions%20Using%20R.pdf
arp=1; %lag length for LP
nstraps=1001; %number of bootstrap replications

[Beff, Sigma] = lp_gls_boot(y,arp,nstraps,maxh,blocksize,badj); 


%Beff is the qxqx(maxh)xnstraps matrix of the wold impulse responses bootstrap replications. For
%example Beff(:,:,h,n)' gives the nth bootstrap draw of the horizon h
%wold impulse (trueirf(:,:,h)) (note that there is a there is a transpose symbol after Beff(:,:,h,n))
%If your interested in all of the LP coefficients and not just the Wold
%impulse responses, add the optional seventh input full where full=1.

%Sigma is the q*q*nstraps matrix of the covariance matrix bootrap
%replications.


%To calculate the 95% conf intervals for the horizon 5 wold
%impulse responses you can do the following
p=[.025 .975];
LPHolder=zeros(q*q,nstraps);
h=5;

%vectorizing the horizon
for n=1:nstraps 
LPHolder(:,n)=reshape(Beff(:,:,h,n)',[],1);
end

       lp_intervals= quantile(LPHolder',p)'; % 95% conf intervals for horizon 5 wold impulse responses

       %The true wold impulse respsonses at horizon 5 are
       reshape(trueirf(:,:,h),[],1);

