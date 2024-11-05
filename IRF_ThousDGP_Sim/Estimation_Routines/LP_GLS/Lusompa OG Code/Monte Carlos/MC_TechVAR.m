function [CoverageTechVAR, LengthTechVAR, BiasTechVAR] = MC_TechVAR(maxh,M,nstraps,long_auto)


%Monte Carlo Code for the Technology VAR in the paper.
%Calculates coverage, average length, and average bias.


%Input:
    %maxh is maximum horizon to be calculated for coverage, average length, and bias
    %M is the number of Monte Carlo simulations to run
    %nstraps is the number of bootstrap replications for the bootstrap estimators
    %long_auto is the lag length to be used when estimating the VAR to be used when generating data (in the paper this is 16 b4 robustness checks)

%Output:
    %CoverageTechVAR is a (q*q)x(maxh) matrix of coverage estimates over all M Monte Carlos
    %LengthTechVAR is a  (q*q)x(maxh) matrix of average length over all M Monte Carlos
    %BiasTechVAR is a (q*q)x(maxh) matrix of average bias over all M Monte Carlos


    load('Technology Shock Data MC')

True_data=[realgdppercapitagrowth realstockpricespercapitagrowth laborproductivitygrowth techshockgrowth]';
[q, T]=size(True_data); 
% names=char('Output Growth', 'SP 500 Groth','Labor Produtivity','Tech Shock'); 
% Ynames=cellstr(names);

tic


burnin=5000;
burnin=burnin+long_auto;


[Sim_Data, trueirf] = sim_data_VAR(True_data,maxh,long_auto,burnin,M,1) ;


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

Temp_data=Sim_Data(:,:,m);
saveY=Temp_data; 



% arp=4; %lag length
%uncomment below to do lag length selection instead
maxlag=long_auto;
[qAIC, qHQC, qBIC, qLBQ, qLLR]= lagselect(Temp_data,maxlag,1,0);
arp=qAIC;

             p=1+arp*q; 
             F=zeros(q*arp,T-arp); 
             for j=1:arp, 
                F((j-1)*q+(1:q),:)=Temp_data(:,(arp-j+1):(T-j)); 
             end
             Temp_data(:,1:arp)=[]; Teffective=size(Temp_data,2); 
             F=[ones(1,Teffective);F];

  
  [~,~,~,resids]= ols(Temp_data',F',0) ;
 
  
for i=1:q
[reject,~]=lbqtest(resids(:,i),'Lags',20,'Alpha',0.1);


while reject==1
    Temp_data=saveY; 


arp=arp+1;

             p=1+arp*q; 
             F=zeros(q*arp,T-arp); 
             for j=1:arp, 
                F((j-1)*q+(1:q),:)=Temp_data(:,(arp-j+1):(T-j)); 
             end
             Temp_data(:,1:arp)=[]; Teffective=size(Temp_data,2); 
             F=[ones(1,Teffective);F];
   
             
             
            [~,~,~,resids]= ols(Temp_data',F',0) ;
  
  

[reject,~]=lbqtest(resids(:,i),'Lags',20,'Alpha',0.1);   
             
end
end

    Temp_data=saveY; 


confint=.95;


blocksize=1;
badj=1;
[Beff, ~] = lp_gls_boot(Temp_data,arp,nstraps,maxh,blocksize,badj);  
[CoverageLP_GLS_boot_temp, LengthLP_GLS_boot_temp, BiasLP_GLS_boot_temp] = lp_gls_boot_summary(Beff,maxh,nstraps,confint,trueirf,q);


blocksize=1;
badj=0;
[Beff2, ~] = lp_gls_boot(Temp_data,arp,nstraps,maxh,blocksize,badj);  
[CoverageLP_GLS_boot_temp2, LengthLP_GLS_boot_temp2, BiasLP_GLS_boot_temp2] = lp_gls_boot_summary(Beff2,maxh,nstraps,confint,trueirf,q);


[Beff3, s_errors3] = lp_gls_analytical(Temp_data,arp,maxh,badj) ;
[CoverageLP_GLS_analytic_temp, LengthLP_GLS_analytic_temp, BiasLP_GLS_analytic_temp] = lp_gls_analytical_summary(Beff3,s_errors3,maxh,confint,trueirf,T,arp,q);


[CoverageLPHARtemp, LengthLPHARtemp, BiasLPHARtemp] = lp_har_summary(Temp_data',arp,maxh,confint,trueirf);


badj=1;
[rhohat, ~] = var_boot(Temp_data,arp,nstraps,badj);
[CoverageVARtemp, LengthVARtemp, BiasVARtemp] = var_boot_summary(rhohat,maxh,nstraps,confint,arp,q,trueirf);


[CoverageVAR_analytictemp, LengthVAR_analytictemp, BiasVAR_analytictemp] = var_analytical_summary(Temp_data,arp,maxh,trueirf,confint);



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


     ['Monte Carlo Iteration ',int2str(m),' of ',int2str(M),' of Tech VAR']


end

 CoverageTechVAR=[mean(CoverageLP_GLS_boot,3) mean(CoverageLP_GLS_boot2,3) mean(CoverageLP_GLS_analytic,3) mean(CoverageLPHAR,3) mean(CoverageVAR_boot,3) mean(CoverageVAR_analytic,3)];
 
 LengthTechVAR=[mean(LengthLP_GLS_boot,3) mean(LengthLP_GLS_boot2,3) mean(LengthLP_GLS_analytic,3) mean(LengthLPHAR,3) mean(LengthVAR_boot,3) mean(LengthVAR_analytic,3)];
 
 BiasTechVAR=[mean(BiasLP_GLS_boot,3) mean(BiasLP_GLS_boot2,3) mean(BiasLP_GLS_analytic,3) mean(BiasLPHAR,3) mean(BiasVAR_boot,3) mean(BiasVAR_analytic,3)];

  
  toc