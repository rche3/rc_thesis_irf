function [CoverageSims, LengthSims,BiasSims] = MC_Sims_SS(maxh,M,nstraps,samplesize)

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
A=[0.0552 0.0746 -0.4402 -0.0909 -1.1461 0.1544 0.4798
0.0034 0.9723 -0.0267 -0.0049 -1.4857 0.0185 -0.0037
-0.0603 -0.0278 0.4806 0.0424 -0.0452 0.0307 0.0875
0.0175 -0.0179 -0.1398 0.7656 -0.0410 0.0077 0.0200
0.0000 -0.0000 -0.0000 -0.0000 0.0000 0.0000 0.0000
0.1417 -0.1717 -1.1291 -0.2051 -0.8958 0.7831 -0.1558
0.0319 0.1411 -0.2543 -0.0600 -1.2136 -0.0153 0.6513];

B=[-1.1461 0.3006
-1.4857 0.0111
-0.0452 0.0096
-0.0410 -0.0105
0.0000 1.0000
-0.8958 0.4691
-1.2136 0.2552];

C=[-0.0000 -0.0000 -0.0000 0.0000 1.0000 0.0000 0.0000
0.0319 0.1411 -0.2543 -0.0600 -1.2136 -0.0153 0.6513];

D=[1.0000 -0.0000
-1.2136 0.2552];

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

Errors=diag([.01 .005])*randn(q,T+burnin);

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

  
  CoverageSims= [mean(CoverageLP_GLS_boot,3) mean(CoverageLP_GLS_boot2,3) mean(CoverageLP_GLS_analytic,3) mean(CoverageLPHAR,3) mean(CoverageVAR_boot,3) mean(CoverageVAR_analytic,3)];
 
 LengthSims=[mean(LengthLP_GLS_boot,3) mean(LengthLP_GLS_boot2,3) mean(LengthLP_GLS_analytic,3) mean(LengthLPHAR,3) mean(LengthVAR_boot,3) mean(LengthVAR_analytic,3)];
 
  BiasSims=[mean(BiasLP_GLS_boot,3) mean(BiasLP_GLS_boot2,3) mean(BiasLP_GLS_analytic,3) mean(BiasLPHAR,3) mean(BiasVAR_boot,3) mean(BiasVAR_analytic,3)];

  
  toc