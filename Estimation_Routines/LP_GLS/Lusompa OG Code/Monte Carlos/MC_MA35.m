function [CoverageMA35, LengthMA35,BiasMA35] = MC_MA35(maxh,M,nstraps,samplesize)

%Monte Carlo Code for the MA(35).
%Calculates coverage, average length, and average bias.


%Input:
    %maxh is maximum horizon to be calculated for coverage, average length, and bias
    %M is the number of Monte Carlo simulations to run
    %nstraps is the number of bootstrap replications for the bootstrap estimators
    %samplesize is desired sample size for the data generating process

%Output:
    %Coveragearma is a (q*q)x(maxh) matrix of coverage estimates over all M Monte Carlos
    %Lengtharma is a  (q*q)x(maxh) matrix of average length over all M Monte Carlos
    %Biasarma is a (q*q)x(maxh) matrix of average bias over all M Monte Carlostic

x_bar=35;
a=.5;
b=12;
c=6;
k=zeros(x_bar,1);

for j=1:x_bar
k(j)=a*exp(-(((j-b)/c).^2));
end

k=k./sum(k);

Sigma=1;


modSim = arima('Constant',0,'AR',0,'MA',k,'Variance',Sigma);
trueirf=impulse(modSim,maxh+1);
trueirf=trueirf(2:end);

CoverageLP_GLS_boot=zeros(maxh,M);
CoverageLP_GLS_boot2=zeros(maxh,M);
CoverageLP_GLS_analytic=zeros(maxh,M);
CoverageLPHAR=zeros(maxh,M);
CoverageVAR_boot=zeros(maxh,M);
CoverageVAR_analytic=zeros(maxh,M);


LengthLP_GLS_boot=zeros(maxh,M);
LengthLP_GLS_boot2=zeros(maxh,M);
LengthLP_GLS_analytic=zeros(maxh,M);
LengthLPHAR=zeros(maxh,M);
LengthVAR_boot=zeros(maxh,M);
LengthVAR_analytic=zeros(maxh,M);


BiasLP_GLS_boot=zeros(maxh,M);
BiasLP_GLS_boot2=zeros(maxh,M);
BiasLP_GLS_analytic=zeros(maxh,M);
BiasLPHAR=zeros(maxh,M);
BiasVAR_boot=zeros(maxh,M);
BiasVAR_analytic=zeros(maxh,M);


 

parfor m=1:M
% m=1
T=samplesize;
burnin=1500;
confint=.95;





y= simulate(modSim,T+burnin);
y=y(burnin+1:end);



Y=y(:,1)'; [q, T]=size(Y); 

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
 
  
[reject,~]=lbqtest(resids,'Lags',20,'Alpha',0.1);


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
  
  

[reject,~]=lbqtest(resids,'Lags',20,'Alpha',0.1);   
             
end
%end of lag length selection


blocksize=1;
badj=1;


[Beff, ~] = lp_gls_boot(y',arp,nstraps,maxh,blocksize,badj)  ;
Beff=squeeze(Beff);
[CoverageLP_GLS_boot_temp, LengthLP_GLS_boot_temp, BiasLP_GLS_boot_temp] =  lp_gls_boot_summary(Beff,maxh,nstraps,confint,trueirf,q)  ;


badj=0;

[Beff2, ~] = lp_gls_boot(y',arp,nstraps,maxh,blocksize,badj)  ;
Beff2=squeeze(Beff2);
[CoverageLP_GLS_boot_temp2, LengthLP_GLS_boot_temp2, BiasLP_GLS_boot_temp2] =  lp_gls_boot_summary(Beff2,maxh,nstraps,confint,trueirf,q)  ;


[Beff, s_errors] = lp_gls_analytical(y',arp,maxh,badj) ;
Beff=squeeze(Beff);
s_errors=squeeze(s_errors);
[CoverageLP_GLS_analytic_temp, LengthLP_GLS_analytic_temp, BiasLP_GLS_analytic_temp] = lp_gls_analytical_summary(Beff,s_errors,maxh,confint,trueirf,T,arp,q); 


[CoverageLPHARtemp, LengthLPHARtemp, BiasLPHARtemp] = lp_har_summary(y,arp,maxh,confint,trueirf);


badj=1;

[rhohat, ~] = var_boot(y',arp,nstraps,badj); 
[CoverageVAR_boottemp, LengthVAR_boottemp, BiasVAR_boottemp] = var_boot_summary(rhohat,maxh,nstraps,confint,arp,q,trueirf) ;


[CoverageVAR_analytictemp, LengthVAR_analytictemp, BiasVAR_analytictemp] = var_analytical_infinite_summary(y',arp,maxh,trueirf,confint);


     
CoverageLP_GLS_boot(:,m)=CoverageLP_GLS_boot_temp;
CoverageLP_GLS_boot2(:,m)=CoverageLP_GLS_boot_temp2;
CoverageLP_GLS_analytic(:,m)=CoverageLP_GLS_analytic_temp;
CoverageLPHAR(:,m)=CoverageLPHARtemp;
CoverageVAR_boot(:,m)=CoverageVAR_boottemp;
CoverageVAR_analytic(:,m)=CoverageVAR_analytictemp;


LengthLP_GLS_boot(:,m)=LengthLP_GLS_boot_temp;
LengthLP_GLS_boot2(:,m)=LengthLP_GLS_boot_temp2;
LengthLP_GLS_analytic(:,m)=LengthLP_GLS_analytic_temp;
LengthLPHAR(:,m)=LengthLPHARtemp;
LengthVAR_boot(:,m)=LengthVAR_boottemp;
LengthVAR_analytic(:,m)=LengthVAR_analytictemp;


BiasLP_GLS_boot(:,m)=BiasLP_GLS_boot_temp;
BiasLP_GLS_boot2(:,m)=BiasLP_GLS_boot_temp2;
BiasLP_GLS_analytic(:,m)=BiasLP_GLS_analytic_temp;
BiasLPHAR(:,m)=BiasLPHARtemp;
BiasVAR_boot(:,m)=BiasVAR_boottemp;
BiasVAR_analytic(:,m)=BiasVAR_analytictemp;


     ['Monte Carlo Iteration ',int2str(m),' of ',int2str(M),' For MA(35)']

end


 CoverageMA35=[mean(CoverageLP_GLS_boot,2) mean(CoverageLP_GLS_boot2,2) mean(CoverageLP_GLS_analytic,2) mean(CoverageLPHAR,2) mean(CoverageVAR_boot,2) mean(CoverageVAR_analytic,2)];
 
 LengthMA35=[mean(LengthLP_GLS_boot,2) mean(LengthLP_GLS_boot2,2) mean(LengthLP_GLS_analytic,2) mean(LengthLPHAR,2) mean(LengthVAR_boot,2) mean(LengthVAR_analytic,2)];
 
  BiasMA35=[mean(BiasLP_GLS_boot,2) mean(BiasLP_GLS_boot2,2) mean(BiasLP_GLS_analytic,2) mean(BiasLPHAR,2) mean(BiasVAR_boot,2) mean(BiasVAR_analytic,2)];

  
  toc

