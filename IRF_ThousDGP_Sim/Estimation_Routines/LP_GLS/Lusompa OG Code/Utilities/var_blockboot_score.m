function [rhohat, Sigma_hold] = var_blockboot_score(y,arp,nstraps,blocksize,badj) 

%Score block bootstrap VAR with bias adjustment option.

%Input:
    %y is the qxT matrix of variables
    %arp is the lag length
    %nstraps is the number of bootstrap replications
    %blocksize is the blocksize of the bootstrap
    %badj equals 1 for bias adjustment and 0 otherwise

%Output
    %rhohat is the pxqxnstraps matrix of the VAR coefficient bootstrap replications. 
    %Sigma_hold is the qxqxnstraps matrix of the covariance matrix bootrap horizon 1 LP (VAR) residuals replications.


Y=y; [q, T]=size(Y); 


             p=1+arp*q; % ---- set up TV-VAR(arp)with zero mean
             F=zeros(q*arp,T-arp); 
             for j=1:arp
                F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-j)); 
             end
             Y(:,1:arp)=[]; Teffective=size(Y,2); 
             F=[ones(1,Teffective);F];
             
             

rhohat=zeros(p,q,nstraps);
epshat=zeros(Teffective,q,nstraps);
epshatsq=zeros(q,q,Teffective,nstraps);
   Sigma_hold=zeros(q,q,nstraps);

    samecov=(F*F'/Teffective)\eye(p);

% estimate VAR coefficients:
  
  [betah,~,~,resids]= ols(Y',F',0) ;

    rhohat(:,:,1) =betah;
  epshat(:,:,1) = resids ;  %population residuals for horizon 0
  Sigma_hold(:,:,1)=resids'*resids/Teffective;
  


  %bias adjustment

  if badj==1
        [A,SIGMA,~,~,~]=olsvarc(y(:,1:end)',arp);
  [bcA]=asybc(A,SIGMA,size(y(:,1:end)',1),arp,size(y(:,1:end)',2));

    rhohat(2:end,:,1)=bcA(1:q,:)'; %using bias adjustment coefficients for population
      Sigma_hold(:,:,1)=SIGMA(1:q,1:q);

for t = arp+1:T 
 epshat(t-arp,:,1)=y(:,t) - rhohat(2:end,:,1)'*reshape(y(:,(t-1):-1:(t-arp)),[],1) ;
end 

epshat(:,:,1) =epshat(:,:,1)-mean(epshat(:,:,1));

  end

%end of bias adjustment for population estimates




%bootstrap portion
        
    for n=2:nstraps
        


        numblock=floor(Teffective/blocksize);
        
        eta=randn(numblock,1);
        

        Hold=repmat(eta,1,blocksize);
        Hold=[reshape(Hold',[],1)];
        temp=repmat(randn,1,(Teffective-(numblock*blocksize)));
        Hold=[Hold; reshape(temp,[],1)];
        
%        epshat(:,:,n)=epshat(:,:,1).*Hold;
       epshat(:,:,n)=resids.*Hold;
%          epshat(:,:,n)=datasample(resids,Teffective);

        for t=1:Teffective
       epshatsq(:,:,t,n-1)=((resids(t,:)'*resids(t,:))-Sigma_hold(:,:,1))*Hold(t);
        end
              
                  Sigma_hold(:,:,n)=Sigma_hold(:,:,1)+sum(epshatsq(:,:,:,n-1),3)/Teffective;
                  
                  d = eig(Sigma_hold(:,:,n));
                  isposdef = all(d > 0);
                  
                                    
                  while isposdef==0
                              eta=randn(numblock,1);

                      Hold=repmat(eta,1,blocksize);
        Hold=[reshape(Hold',[],1)];
        temp=repmat(randn,1,(Teffective-(numblock*blocksize)));
        Hold=[Hold; reshape(temp,[],1)];
        
%        epshat(:,:,n)=epshat(:,:,1).*Hold;
       epshat(:,:,n)=resids.*Hold;
%          epshat(:,:,n)=datasample(resids,Teffective);

        for t=1:Teffective
       epshatsq(:,:,t,n-1)=((resids(t,:)'*resids(t,:))-Sigma_hold(:,:,1))*Hold(t);
        end
                      Sigma_hold(:,:,n)=Sigma_hold(:,:,1)+sum(epshatsq(:,:,:,n-1),3)/Teffective;
                  
                  d = eig(Sigma_hold(:,:,n));
                  isposdef = all(d > 0);
                      
                  end
                  
               rhohat(:,:,n)=rhohat(:,:,1)+ samecov*F*epshat(:,:,n)/Teffective ;        

    end

return