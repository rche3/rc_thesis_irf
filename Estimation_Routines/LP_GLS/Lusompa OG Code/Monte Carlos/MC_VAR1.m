function [CoverageVAR1, LengthVAR1,BiasVAR1] = MC_VAR1(maxh,coef,M,nstraps,samplesize)


%Monte Carlo Code for the Kilian and Kim (2011) VAR(1).
%Calculates coverage, average length, and average bias.


%Input:
    %maxh is maximum horizon to be calculated for coverage, average length, and bias
    %coef is the chosen coefficient for the VAR(1) (.97, .9, or .5 in the paper)
    %M is the number of Monte Carlo simulations to run
    %nstraps is the number of bootstrap replications for the bootstrap estimators
    %samplesize is desired sample size for the data generating process

%Output:
    %CoverageVAR1 is a (q*q)x(maxh) matrix of coverage estimates over all M Monte Carlos
    %LengthVAR1 is a  (q*q)x(maxh) matrix of average length over all M Monte Carlos
    %BiasVAR1 is a (q*q)x(maxh) matrix of average bias over all M Monte Carlos

tic


B1=[coef 0;.5 .5];


q=size(B1,1);
trueirf=zeros(q,q,maxh);

for h=1:maxh
    trueirf(:,:,h)=(B1^h);
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
    burnin=1000;


Errors=randn(q,T+burnin);
y=zeros(q,T+burnin);
y(:,1)=randn(q,1);


for t=2:T+burnin
    y(:,t)=B1*y(:,t-1)+Errors(:,t);
end
y=y(:,burnin+1:end);

Y=y; [q, T]=size(Y); 


saveY=Y; 



%lag length selection
maxlag=8;
[qAIC, qHQC, qBIC, qLBQ, qLLR]= lagselect(Y,maxlag,1,0);
arp=qAIC;
laglength(m)=arp;

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


[Beff3, s_errors3] = lp_gls_analytical(y,arp,maxh,badj) ;
[CoverageLP_GLS_analytic_temp, LengthLP_GLS_analytic_temp, BiasLP_GLS_analytic_temp] = lp_gls_analytical_summary(Beff3,s_errors3,maxh,confint,trueirf,T,arp,q);


[CoverageLPHARtemp, LengthLPHARtemp, BiasLPHARtemp] = lp_har_summary(y',arp,maxh,confint,trueirf);


badj=1;
[rhohat, ~] = var_boot(y,arp,nstraps,badj);
[CoverageVARtemp, LengthVARtemp, BiasVARtemp] = var_boot_summary(rhohat,maxh,nstraps,confint,arp,q,trueirf);

[CoverageVAR_analytictemp, LengthVAR_analytictemp, BiasVAR_analytictemp] = var_analytical_summary(y,arp,maxh,trueirf,confint);



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


     ['Monte Carlo Iteration ',int2str(m),' of ',int2str(M),' For VAR(1) a_11 = ',num2str(coef)]

end

 
 CoverageVAR1=[mean(CoverageLP_GLS_boot,3) mean(CoverageLP_GLS_boot2,3) mean(CoverageLP_GLS_analytic,3) mean(CoverageLPHAR,3) mean(CoverageVAR_boot,3) mean(CoverageVAR_analytic,3)];
 
 LengthVAR1=[mean(LengthLP_GLS_boot,3) mean(LengthLP_GLS_boot2,3) mean(LengthLP_GLS_analytic,3) mean(LengthLPHAR,3) mean(LengthVAR_boot,3) mean(LengthVAR_analytic,3)];
 
  BiasVAR1=[mean(BiasLP_GLS_boot,3) mean(BiasLP_GLS_boot2,3) mean(BiasLP_GLS_analytic,3) mean(BiasLPHAR,3) mean(BiasVAR_boot,3) mean(BiasVAR_analytic,3)];

  
  toc

  