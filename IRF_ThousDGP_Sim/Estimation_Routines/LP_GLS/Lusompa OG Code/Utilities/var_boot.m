function [rhohat, Sigma_hold] = var_boot(y,arp,nstraps,badj) 

%VAR Bootstrap Code (recursive design iid)

%Input:
    %y is the qxT matrix of variables
    %arp is the lag length
    %nstraps is the number of bootstrap replications
    %badj equals 1 for bias adjustment and 0 otherwise

%Output
    %rhohat is the pxqxnstraps matrix of the VAR coefficient bootstrap replications. 
    %Sigma_hold is the qxqxnstraps matrix of the covariance matrix bootrap horizon 1 LP (VAR) residuals replications.


Y=y; [q, T]=size(Y); 

y=cat(3,y,zeros(q,T,nstraps-1));


             p=1+arp*q; 
             F=zeros(q*arp,T-arp); 
             for j=1:arp
                F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-j)); 
             end
             Y(:,1:arp)=[]; Teffective=size(Y,2); 
             F=[ones(1,Teffective);F];
             
             

rhohat=zeros(p,q,nstraps);
epshat=zeros(Teffective,q,nstraps);  
   Sigma_hold=zeros(q,q,nstraps);


% estimate VAR coefficients:
  
  [betah,~,~,resids]= ols(Y',F',0) ;

    rhohat(:,:,1) =betah;
  epshat(:,:,1) = resids ;  %population residuals for horizon 0
  Sigma_hold(:,:,1)=resids'*resids/Teffective;
  


  %bias adjustment

  if badj==1
        [A,SIGMA,~,~,~]=olsvarc(y(:,1:end,1)',arp);
  [bcA]=asybc(A,SIGMA,size(y(:,1:end,1)',1),arp,size(y(:,1:end,1)',2));

    rhohat(2:end,:,1)=bcA(1:q,:)'; %using bias adjustment coefficients for population
      Sigma_hold(:,:,1)=SIGMA(1:q,1:q);

      for t = arp+1:T 
 epshat(t-arp,:,1)=y(:,t,1) - rhohat(2:end,:,1)'*reshape(y(:,(t-1):-1:(t-arp),1),[],1) ;
end 

epshat(:,:,1) =epshat(:,:,1)-mean(epshat(:,:,1));

  end

%end of bias adjustment for population estimates




%bootstrap portion
        
                        randomat=randi(T-arp+1,nstraps-1,1);   

for n = 2:nstraps 
    
                 rando=randomat(n-1);
y(:,1:arp,n)=y(:,rando:rando+arp-1,1);  %using random block as initial obs
    
    
    
YourVector=datasample(resids,Teffective);
YourVector=YourVector-mean(YourVector); 
    
for t = arp+1:T 
 y(:,t,n) = rhohat(1,:,1)'+rhohat(2:end,:,1)'*reshape(y(:,(t-1):-1:(t-arp),n),[],1) +YourVector(t-arp,:)' ;
end 
end


for n = 2:nstraps 
    Y=y(:,:,n); [q, T]=size(Y); 

             F=zeros(q*arp,T-arp); 
             for j=1:arp 
                F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-j)); 
             end
             Y(:,1:arp)=[]; Teffective=size(Y,2); 
             F=[ones(1,Teffective);F];
       
  [betah,~,~,resids]= ols(Y',F',0) ;
  rhohat(:,:,n) = betah  ;
  epshat(:,:,n) = resids ;
  epshat(:,:,n) =epshat(:,:,n)-mean(epshat(:,:,n));

    Sigma_hold(:,:,n)=epshat(:,:,n)'*epshat(:,:,n)/Teffective;

%bias adjustment of bootstrap draws (bias adjustment after bias adjustment discussed in Kilian (1998))  
  if badj==1
    [A,SIGMA,~,~,~]=olsvarc(y(:,1:end,n)',arp);
  [bcA]=asybc(A,SIGMA,size(y(:,1:end,n)',1),arp,size(y(:,1:end,n)',2));
  
  rhohat(2:end,:,n)=bcA(1:q,:)';

  for t = arp+1:T 
 epshat(t-arp,:,n)=y(:,t,n) - rhohat(2:end,:,n)'*reshape(y(:,(t-1):-1:(t-arp),n),[],1) ;  %%?
end 

epshat(:,:,n) =epshat(:,:,n)-mean(epshat(:,:,n));

Sigma_hold(:,:,n)=epshat(:,:,n)'*epshat(:,:,n)/Teffective;
  end


end 

return