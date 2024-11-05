function [Sim_Data, trueirf] = sim_data_VAR2(y,maxh,long_auto,burnin,M,badj,samplesize) 

%Estimates a VAR (with bias adjustment if desired) and then generates data
%along with the true irfs (where the original estimate VAR coefficients are treated
%as the true coefficient)

%Input:
    %y is the qxT matrix of variables
    %maxh is the max horizon to be calculated for the trueirf
    %long_auto is the lag length of the estimated VAR
    %burnin is the burnin sample size (doesn't assume initial sample is from stationary distribution)    
    %M is the number of generated samples
    %badj equals 1 for bias adjustment and 0 otherwise
    %samplesize the desired sample size for the simulated data

%Output:
    %Sim_Data is the qxsamplesizexM matrix of simulated data
    %trueirf is the qxqx(maxh) matrix the Wold Coefficients (treating the estimated VAR as the true data generating process)


Y=y; [q, T]=size(Y); 


arp=long_auto;

             p=1+arp*q; 
             F=zeros(q*arp,T-arp); 
             for j=1:arp
                F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-j)); 
             end
             Y(:,1:arp)=[]; Teffective=size(Y,2); 
             F=[ones(1,Teffective);F]; %Design matrix

% estimate VAR coefficients:
  
  [betah,~,~,resids]= ols(Y',F',0) ;

  rhohat=zeros(p,q,1);
  rhohat(:,:,1) =betah;

    Sigma_hold=resids'*resids/Teffective;


    if badj==1 %bias adjustment
  [A,SIGMA,~,~,~]=olsvarc(y(:,1:end,1)',arp);

  [bcA]=asybc(A,SIGMA,size(y(:,1:end,1)',1),arp,size(y(:,1:end,1)',2));

    rhohat(2:end,:,1)=bcA(1:q,:)';
    end


  %simulating data      
Sim_Data=zeros(q,samplesize+burnin,M);

CholSigma=chol(Sigma_hold,'lower');
    
   for m=1:M
for t = arp+1:samplesize+burnin 
 Sim_Data(:,t,m) = rhohat(1,:,1)'+rhohat(2:end,:,1)'*reshape(Sim_Data(:,(t-1):-1:(t-arp),m),[],1) +CholSigma*randn(q,1) ;
end 
   end

   %removing burnin
Sim_Data=Sim_Data(:,burnin+1:end,:);


   %calculating Wold coefficients
trueirf=zeros(q,q,maxh);

Ap=rhohat(2:end,:,1)';
VAR_Companion=[Ap; eye((arp-1)*q) zeros(((arp-1)*q),q)];

for h=1:maxh    
    Hold=VAR_Companion^h;
trueirf(:,:,h)=Hold(1:q,1:q);
end



return