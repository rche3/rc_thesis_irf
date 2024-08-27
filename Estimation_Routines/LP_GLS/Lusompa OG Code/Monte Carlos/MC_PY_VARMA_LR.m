function [CoveragePYLR, LengthPYLR, BiasPYLR] = MC_PY_VARMA_LR(maxh,maxhLR,M,nstraps,samplesize)

%Monte Carlo Code for the Poskit and Yao (2017) VARMA(1,1) with long run restrictions imposed.
%Calculates coverage, average length, and average bias.


%Input:
    %maxh is maximum horizon to be calculated for coverage, average length, and bias
    %M is the number of Monte Carlo simulations to run
    %nstraps is the number of bootstrap replications for the bootstrap estimators
    %samplesize is desired sample size for data generating process

%Output:
    %CoverageVARMA is a (q*q)x(maxh) matrix of coverage estimates over all M Monte Carlos
    %LengthVARMA is a (q*q)x(maxh) matrix of average length over all M Monte Carlos
    %BiasVARMA is a (q*q)x(maxh) matrix of average bias over all M Monte Carlos

    tic


B1=[0.9413 1.0446;0.00060 0.8045];
M1=[-0.2498 -0.9173;-0.1924 -0.7065];

Phi1=eye(2)-B1;
Theta1=eye(2)+M1;

Gamma1=Phi1\Theta1;


Sigma=[0.51863008 0.40575031;0.40575031 0.40089860]*(10^(-3));

GG=Gamma1*Sigma*Gamma1';

P=chol(Sigma,'lower');
G=chol(GG,'lower'); 
Lambda=Gamma1\G;

q=size(B1,1);
trueirf=zeros(q,q,maxh);

for h=1:maxh
    trueirf(:,:,h)=((B1^h)+((B1^(h-1))*M1))*Lambda;
end


CoverageLP_GLS_boot=zeros(maxh,q*q,M);
CoverageLP_GLS_boot2=zeros(maxh,q*q,M);
CoverageLP_GLS_analytic=zeros(maxh,q*q,M);
CoverageLPHAR=zeros(maxh,q*q,M);
CoverageVAR_boot=zeros(maxh,q*q,M);
CoverageVAR_analytic=zeros(maxh,q*q,M);

    
LengthLP_GLS_boot=zeros(maxh,q*q,M);
LengthLP_GLS_boot2=zeros(maxh,q*q,M);
LengthLP_GLS_analytic=zeros(maxh,q*q,M);
LengthLPHAR=zeros(maxh,q*q,M);
LengthVAR_boot=zeros(maxh,q*q,M);
LengthVAR_analytic=zeros(maxh,q*q,M);


BiasLP_GLS_boot=zeros(maxh,q*q,M);
BiasLP_GLS_boot2=zeros(maxh,q*q,M);
BiasLP_GLS_analytic=zeros(maxh,q*q,M);
BiasLPHAR=zeros(maxh,q*q,M);
BiasVAR_boot=zeros(maxh,q*q,M);
BiasVAR_analytic=zeros(maxh,q*q,M);

 
parfor m=1:M
% m=1
q=size(B1,1);


    T=samplesize;
    burnin=5000;


Errors=P*randn(q,T+burnin);
y=zeros(q,T+burnin);
y(:,1)=P*randn(q,1)+M1*P*randn(q,1);


for t=2:T+burnin
    y(:,t)=B1*y(:,t-1)+Errors(:,t)+M1*Errors(:,t-1);
end
y=y(:,burnin+1:end);


arp=8;%lag length; comment out below for lag length selection
% Y=y; [q T]=size(Y); 
% 
% 
% saveY=Y; 
% 
% 
% maxlag=8;
% [qAIC, qHQC, qBIC, qLBQ, qLLR]= lagselect(Y,maxlag,1,0);
% %could need to do lag length selection around here
% arp=qAIC;
% 
%              p=1+arp*q; 
%              F=zeros(q*arp,T-arp); 
%              for j=1:arp, 
%                 F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-j)); 
%              end
%              Y(:,1:arp)=[]; Teffective=size(Y,2); 
%              F=[ones(1,Teffective);F];
%              

%   
%   [~,~,~,resids]= ols(Y',F',0) ;
%  
%   
% for i=1:q
% [reject,~]=lbqtest(resids(:,i),'Lags',20,'Alpha',0.1);
% 
% 
% while reject==1
%     Y=saveY; 
% 
% 
% arp=arp+1;
% 
%              p=1+arp*q; 
%              F=zeros(q*arp,T-arp); 
%              for j=1:arp, 
%                 F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-j)); 
%              end
%              Y(:,1:arp)=[]; Teffective=size(Y,2); 
%              F=[ones(1,Teffective);F];
%    
%              
%              
%             [~,~,~,resids]= ols(Y',F',0) ;
%   
%   
% 
% [reject,~]=lbqtest(resids(:,i),'Lags',20,'Alpha',0.1);   
%              
% end
% end
%%end of lag length selection



confint=.95;



blocksize=1;
badj=1;


[Beff, Sigma_hold] = lp_gls_boot(y,arp,nstraps,maxhLR,blocksize,badj);  
[CoverageLP_GLS_boot_temp, LengthLP_GLS_boot_temp, BiasLP_GLS_boot_temp] = lp_gls_boot_summary_longrun(Beff,Sigma_hold,maxh,maxhLR,nstraps,confint,trueirf,q);


badj=0;
[Beff2, Sigma_hold2] = lp_gls_boot(y,arp,nstraps,maxhLR,blocksize,badj);  
[CoverageLP_GLS_boot_temp2, LengthLP_GLS_boot_temp2, BiasLP_GLS_boot_temp2] = lp_gls_boot_summary_longrun(Beff2,Sigma_hold2,maxh,maxhLR,nstraps,confint,trueirf,q);



[Beff3, Bigcov,Sigma_hold,epshat] = lp_gls_analytical_longrun(y,arp,maxhLR,badj) ;
[CoverageLP_GLS_analytic_temp, LengthLP_GLS_analytic_temp, BiasLP_GLS_analytic_temp] = lp_gls_analytical_summary_longrun(Beff3,Bigcov,Sigma_hold,epshat,maxh,maxhLR,confint,trueirf,T,arp,q);


[CoverageLPHARtemp, LengthLPHARtemp, BiasLPHARtemp] = lp_har_summary_longrun(y,arp,maxh,maxhLR,confint,trueirf);


badj=1;
[rhohat, Sigma_holdrho] = var_boot(y,arp,nstraps,badj);
[CoverageVARtemp, LengthVARtemp, BiasVARtemp] = var_boot_summary_longrun(rhohat,Sigma_holdrho,maxh,nstraps,confint,arp,q,trueirf);


[CoverageVAR_analytictemp, LengthVAR_analytictemp, BiasVAR_analytictemp] = var_analytical_infinite_summary_longrun(y,arp,maxh,trueirf,confint);



CoverageLP_GLS_boot(:,:,m)=CoverageLP_GLS_boot_temp;
CoverageLP_GLS_boot2(:,:,m)=CoverageLP_GLS_boot_temp2;
CoverageLP_GLS_analytic(:,:,m)=CoverageLP_GLS_analytic_temp;
CoverageLPHAR(:,:,m)=CoverageLPHARtemp;
CoverageVAR_boot(:,:,m)=CoverageVARtemp;
CoverageVAR_analytic(:,:,m)=CoverageVAR_analytictemp;


LengthLP_GLS_boot(:,:,m)=LengthLP_GLS_boot_temp;
LengthLP_GLS_boot2(:,:,m)=LengthLP_GLS_boot_temp2;
LengthLP_GLS_analytic(:,:,m)=LengthLP_GLS_analytic_temp;
LengthLPHAR(:,:,m)=LengthLPHARtemp;
LengthVAR_boot(:,:,m)=LengthVARtemp;
LengthVAR_analytic(:,:,m)=LengthVAR_analytictemp;


BiasLP_GLS_boot(:,:,m)=BiasLP_GLS_boot_temp;
BiasLP_GLS_boot2(:,:,m)=BiasLP_GLS_boot_temp2;
BiasLP_GLS_analytic(:,:,m)=BiasLP_GLS_analytic_temp;
BiasLPHAR(:,:,m)=BiasLPHARtemp;
BiasVAR_boot(:,:,m)=BiasVARtemp;
BiasVAR_analytic(:,:,m)=BiasVAR_analytictemp;


     ['Monte Carlo Iteration ',int2str(m),' of ',int2str(M),' for VARMA(1,1) PK LR restrictions imposed']


end

  
CoveragePYLR=[mean(CoverageLP_GLS_boot,3) mean(CoverageLP_GLS_boot2,3) mean(CoverageLP_GLS_analytic,3) mean(CoverageLPHAR,3) mean(CoverageVAR_boot,3) mean(CoverageVAR_analytic,3)];
 
LengthPYLR=[mean(LengthLP_GLS_boot,3) mean(LengthLP_GLS_boot2,3) mean(LengthLP_GLS_analytic,3) mean(LengthLPHAR,3) mean(LengthVAR_boot,3) mean(LengthVAR_analytic,3)];
 
BiasPYLR=[mean(BiasLP_GLS_boot,3) mean(BiasLP_GLS_boot2,3) mean(BiasLP_GLS_analytic,3) mean(BiasLPHAR,3) mean(BiasVAR_boot,3) mean(BiasVAR_analytic,3)];

  
  toc