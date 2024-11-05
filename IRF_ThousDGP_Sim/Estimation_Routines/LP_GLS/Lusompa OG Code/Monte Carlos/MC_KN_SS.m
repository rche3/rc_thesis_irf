function [CoverageKN, LengthKN,BiasKN] = MC_KN_SS(maxh,M,nstraps,samplesize)

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

A=[.9 0 0;0 .95 0;.545 0 .5143];

B=[1 0 0;0 1 0;.6055 0 .6858];

C=[.545 0 .5143; 1.3377 .95 -.8258;1.3418 0 -.5596];

D=[.6055 0 .6858;1.4863 1 -1.1011;1.4909 0 -.7462];

mid=A-(B/D)*C;
back=B/D;



q=size(C,1);
varlagcoeff=zeros(q,q,maxh);
trueirf=zeros(q,q,maxh);

    varlagcoeff(:,:,1)=C*back;
    trueirf(:,:,1)=varlagcoeff(:,:,1);

for h=2:maxh

        varlagcoeff(:,:,h)=C*(mid^(h-1))*back;

    for j=1:h-1
    trueirf(:,:,h)=trueirf(:,:,h)+trueirf(:,:,h-j)*varlagcoeff(:,:,j);
    end

    trueirf(:,:,h)=trueirf(:,:,h)+varlagcoeff(:,:,h);
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

q=size(C,1);

    T=samplesize;
    burnin=5000;

Errors=diag([.3/100 .6/100 .2/100])*randn(q,T+burnin);

y=zeros(q,T+burnin);
s=zeros(size(A,1),T+burnin);

s(:,1)=B*Errors(:,1);

for t=2:T+burnin
    s(:,t)=A*s(:,t-1)+B*Errors(:,t);
        y(:,t)=C*s(:,t-1)+D*Errors(:,t);
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


     ['Monte Carlo Iteration ',int2str(m),' of ',int2str(M),' for Sims SS Model']


end

  
  CoverageKN= [mean(CoverageLP_GLS_boot,3) mean(CoverageLP_GLS_boot2,3) mean(CoverageLP_GLS_analytic,3) mean(CoverageLPHAR,3) mean(CoverageVAR_boot,3) mean(CoverageVAR_analytic,3)];
 
 LengthKN=[mean(LengthLP_GLS_boot,3) mean(LengthLP_GLS_boot2,3) mean(LengthLP_GLS_analytic,3) mean(LengthLPHAR,3) mean(LengthVAR_boot,3) mean(LengthVAR_analytic,3)];
 
  BiasKN=[mean(BiasLP_GLS_boot,3) mean(BiasLP_GLS_boot2,3) mean(BiasLP_GLS_analytic,3) mean(BiasLPHAR,3) mean(BiasVAR_boot,3) mean(BiasVAR_analytic,3)];

  
  toc