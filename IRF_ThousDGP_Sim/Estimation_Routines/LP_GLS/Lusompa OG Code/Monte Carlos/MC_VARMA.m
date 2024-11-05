function [CoverageVARMA, LengthVARMA,BiasVARMA] = MC_VARMA(maxh,M,nstraps,samplesize)

%Monte Carlo Code for the Kilian and Kim (2011) VARMA(1,1).
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


B1=[0.5417 -0.1971 -0.9395
    0.0400   0.9677   0.0323
    -0.0015  0.0829   0.8080];

M1=[-0.1428 -1.5133 -0.705
    -0.0202  0.0309   0.1561
    0.0227   0.1178 -0.0153];

P=[9.235200 0 0
   -1.4343  3.607 0
   -0.7756  1.2296  2.7555];

q=size(B1,1);
trueirf=zeros(q,q,maxh);

for h=1:maxh
    trueirf(:,:,h)=(B1^h)+((B1^(h-1))*M1);
end



CoverageLP_GLS_boot=zeros(maxh,q*q,M);
CoverageLP_GLS_boot2=zeros(maxh,q*q,M);
CoverageLP_GLS_analytic=zeros(maxh,q*q,M);
CoverageLPHAR=zeros(maxh,q*q,M);
CoverageVAR_boot=zeros(maxh,q*q,M);

LengthLP_GLS_boot=zeros(maxh,q*q,M);
LengthLP_GLS_boot2=zeros(maxh,q*q,M);
LengthLP_GLS_analytic=zeros(maxh,q*q,M);
LengthLPHAR=zeros(maxh,q*q,M);
LengthVAR_boot=zeros(maxh,q*q,M);

BiasLP_GLS_boot=zeros(maxh,q*q,M);
BiasLP_GLS_boot2=zeros(maxh,q*q,M);
BiasLP_GLS_analytic=zeros(maxh,q*q,M);
BiasLPHAR=zeros(maxh,q*q,M);
BiasVAR_boot=zeros(maxh,q*q,M);

 
parfor m=1:M
% m=1
q=size(B1,1);


    T=samplesize;
    burnin=1000;


Errors=P*randn(q,T+burnin);
y=zeros(q,T+burnin);
y(:,1)=P*randn(q,1)+M1*P*randn(q,1);


for t=2:T+burnin
    y(:,t)=B1*y(:,t-1)+Errors(:,t)+M1*Errors(:,t-1);
end
y=y(:,burnin+1:end);

Y=y; q=size(Y,1); 


saveY=Y; 



%lag length selection
maxlag=8;
[qAIC, qHQC, qBIC, qLBQ, qLLR]= lagselect(Y,maxlag,1,0);
arp=qAIC;
             p=1+arp*q; 
             F=zeros(q*arp,T-arp); 
             for j=1:arp
                F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-j)); 
             end
             Y(:,1:arp)=[]; Teffective=size(Y,2); 
             F=[ones(1,Teffective);F];
             
  
  [~,~,~,resids]= ols(Y',F',0) ;
 
  
for i=1:q
[reject,~]=lbqtest(resids(:,i),'Lags',20,'Alpha',0.1);


while reject==1
    Y=saveY; 


arp=arp+1;

             p=1+arp*q;
             F=zeros(q*arp,T-arp); 
             for j=1:arp
                F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-j)); 
             end
             Y(:,1:arp)=[]; Teffective=size(Y,2); 
             F=[ones(1,Teffective);F];
   
             
             
            [~,~,~,resids]= ols(Y',F',0) ;
  
  

[reject,~]=lbqtest(resids(:,i),'Lags',20,'Alpha',0.1);   
             
end
end
%end of lag length selection




confint=.95;



blocksize=1;
badj=1;

[Beff, ~] = lp_gls_boot(y,arp,nstraps,maxh,blocksize,badj);  
[CoverageLP_GLS_boot_temp, LengthLP_GLS_boot_temp, BiasLP_GLS_boot_temp] = lp_gls_boot_summary(Beff,maxh,nstraps,confint,trueirf,q);


badj=0;


[Beff2, ~] = lp_gls_boot(y,arp,nstraps,maxh,blocksize,badj);  

[CoverageLP_GLS_boot_temp2, LengthLP_GLS_boot_temp2, BiasLP_GLS_boot_temp2] = lp_gls_boot_summary(Beff2,maxh,nstraps,confint,trueirf,q);


[Beff2, s_errors2] = lp_gls_analytical(y,arp,maxh,badj) ;
[CoverageLP_GLS_analytic_temp, LengthLP_GLS_analytic_temp, BiasLP_GLS_analytic_temp] = lp_gls_analytical_summary(Beff2,s_errors2,maxh,confint,trueirf,T,arp,q);

[CoverageLPHARtemp, LengthLPHARtemp, BiasLPHARtemp] = lp_har_summary(y',arp,maxh,confint,trueirf);

badj=1;

[rhohat, ~] = var_boot(y,arp,nstraps,badj);
[CoverageVARtemp, LengthVARtemp, BiasVARtemp] = var_boot_summary(rhohat,maxh,nstraps,confint,arp,q,trueirf);

[CoverageVAR_analytictemp, LengthVAR_analytictemp, BiasVAR_analytictemp] = var_analytical_infinite_summary(y,arp,maxh,trueirf,confint);



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


     ['Monte Carlo Iteration ',int2str(m),' of ',int2str(M),' for VARMA(1,1) Kilian and Kim']


end

  
  CoverageVARMA= [mean(CoverageLP_GLS_boot,3) mean(CoverageLP_GLS_boot2,3) mean(CoverageLP_GLS_analytic,3) mean(CoverageLPHAR,3) mean(CoverageVAR_boot,3) mean(CoverageVAR_analytic,3)];
 
 LengthVARMA=[mean(LengthLP_GLS_boot,3) mean(LengthLP_GLS_boot2,3) mean(LengthLP_GLS_analytic,3) mean(LengthLPHAR,3) mean(LengthVAR_boot,3) mean(LengthVAR_analytic,3)];
 
  BiasVARMA=[mean(BiasLP_GLS_boot,3) mean(BiasLP_GLS_boot2,3) mean(BiasLP_GLS_analytic,3) mean(BiasLPHAR,3) mean(BiasVAR_boot,3) mean(BiasVAR_analytic,3)];

  
  toc