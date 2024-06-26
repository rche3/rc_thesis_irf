function [Beff, Sigma_hold] = lusompa_lp_gls_boot(y,arp,nstraps,maxh,blocksize,badj,full) 

%LP GLS Bootstrap Code
cd /Users/rogerchen/Documents/MATLAB/Thesis_24/LPW_24_v3/lp_var_simul/Estimation_Routines/LP_GLS/
%Input:
    %y is the qxT matrix of variables
    %arp is the lag length
    %nstraps is the number of bootstrap replications
    %maxh is the maximum desired LP horizon
    %blocksize is the blocksize of the bootstrap
    %badj equals 1 for bias adjustment and 0 otherwise

%Optional Input:
    %full is an optional input argument. If full is 0 it gives Wold
        %coefficients and if full=1 it gives all of the LP coefficients
        %(not just the wold). If full is not included in the input argument
        %it'll give you just the Wold coefficients.

%Output:
    %Beff is the qxqx(maxh)xnstraps matrix of the wold impulse responses bootstrap replications (default). For
            %example Beff(:,:,h,n)' gives the nth bootstrap draw of the
            %horizon h Wold impulse. If full is 1, then it'll give you all
            %of the LP coefficients (as specified in the paper) not just the Wold estimates.
    %Sigma_hold is the qxqxnstraps matrix of the covariance matrix bootrap horizon 1 LP (VAR) residuals replications.

    if nargin<7
        full=0;
    end

Y=y; [q, T]=size(Y); 
  
 
%population estimate for h=1
for h = 1 
      
saveY=Y; 
  
    p=1+arp*q; 
F=zeros(q*arp,T-arp-h+1); 
for j=1:arp 
    F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-h+1-j));  
end
  
saveF=F;
  
yeff=saveY(:,arp+h:end)'; %Dependent variable
 
Teffective=size(yeff',2);
F=[ones(1,Teffective);saveF]; % Design Matrix

[betahat,~,~,resids] = ols(yeff,F',0) ;
 
 
epshat=zeros(Teffective,q,nstraps); %matrix for bootstrap residuals
epshatsq=zeros(q,q,Teffective,nstraps-1); 


epshat(:,:,1)=resids;

    Beff = zeros(p,q,maxh,nstraps) ; %matrix for LP bootstrap replications
    Sigma_hold=zeros(q,q,nstraps); %matrix for bootstrap covariance replications
%     Autohold=zeros(p,q,2*T,h,nstraps-1);

 
    Beff(:,:,h,1) = betahat ;
            Sigma_hold(:,:,1)=resids'*resids/Teffective;
    
    
        saveF_h1=F;
    samecov=(saveF_h1*saveF_h1'/Teffective)\eye(p); %Gamma inverse in the paper
    
    Qsome=zeros(p,p,maxh-1);
 
    for l=1:maxh-1
        Qsome(:,:,l)=(saveF_h1(:,1+l:end)*saveF_h1(:,1:end-l)'/(Teffective)); %autocovariances of X in the paper
    end

 




        %bias adjustment for h=1
if badj==1

    biaslp=zeros(p-1,q);

[A,SIGMA,~,~,~]=olsvarc(y(:,1:end-h)',arp);

phitwid=A;
omegaUtwid=SIGMA;
nZtwid=size(SIGMA,1);
nk=nZtwid;
PX=eye(nZtwid);

%calculating bias for each equation
for dim=1:q

EetaZtwid0=(F(2:end,arp+1:end)*resids(1:end-arp,dim))'/(Teffective-arp);

[vbias] = proc_vb_ma0(nZtwid, phitwid, omegaUtwid, nk, PX, EetaZtwid0); %bias adjustment function

biaslp(:,dim)=vbias;

end


biaslp= biaslp./Teffective; %estimate of bias
betahat_adj=betahat(2:end,:)-biaslp;
 
 
Apbefore=betahat(2:end,:)';
VAR_Companionbefore=[Apbefore; eye((arp-1)*q) zeros(((arp-1)*q),q)];
 
 
%restricting bias to stationarity
if any(abs(eig(VAR_Companionbefore))>=1)
    Beff(2:end,:,h,1) = betahat(2:end,:) ;
else
    Beff(2:end,:,h,1) = betahat_adj ;
end
 
Ap=Beff(2:end,:,h,1)';
VAR_Companion=[Ap; eye((arp-1)*q) zeros(((arp-1)*q),q)];
 
john=0;
 
while all(abs(eig(VAR_Companionbefore))<1) && any(abs(eig(VAR_Companion))>=1)
john=john+1;
    Beff(2:end,:,h,1)=betahat(2:end,:)-(biaslp*(1-.01*john));
    
if john==100
        break
end  
    
    Ap=Beff(2:end,:,h,1)';
VAR_Companion=[Ap; eye((arp-1)*q) zeros(((arp-1)*q),q)];

 
end
%%%%% end of bias adjustment
 

      epshat(:,:,1) = yeff-F'*Beff(:,:,1,1);  %population residuals for horizon 1      
      epshat(:,:,1)=epshat(:,:,1)-mean(epshat(:,:,1));
      
end

 
  %generate bootstrap samples for horizon 1 LP
  
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
      
          Beff(:,:,h,n)=Beff(:,:,h,1)+ samecov*saveF_h1*epshat(:,:,n)/Teffective ;

    end
 
end
 
 
%population estimate for h>=2
  for h = 2:maxh 
      
          epscorrect = zeros(T-arp-h+1,q) ; 
          
    %gls correction for population estimate
    for hsub = 1:h-1 
      epscorrect = epscorrect +  epshat(h-hsub:T-arp-hsub,:,1)*Beff(2:q+1,:,hsub,1) ;
    end 
        
    Y=y;  
saveY=Y; 
  
yeff=saveY(:,arp+h:end)'-epscorrect; %Dependent Variable


    p=1+arp*q; 
F=zeros(q*arp,T-arp-h+1); 
for j=1:arp 
    F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-h+1-j)); 
 
end
 
saveF=F; 
 
 
 
Teffective=size(yeff',2);
F=[ones(1,Teffective);saveF]; %Desing Matrix
 
 
    %population estimate
    [betahat,~,~,resids] = ols(yeff,F',0) ;
    
    Beff(:,:,h,1) = betahat ;
 
    
 
    %bias correction
[A,SIGMA,~,~,~]=olsvarc(y(:,1:end-h)',arp);

phitwid=A;
omegaUtwid=SIGMA;
nZtwid=size(SIGMA,1);
nk=nZtwid;
PX=eye(nZtwid);



biaslp=zeros(p-1,q);

if badj==1

for dim=1:q

EetaZtwid0=(F(2:end,arp+1:end)*resids(1:end-arp,dim))'/(Teffective-arp);

[vbias] = proc_vb_ma0(nZtwid, phitwid, omegaUtwid, nk, PX, EetaZtwid0);

biaslp(:,dim)=vbias;

end
  

%%%%%%%%%%%% restricted the bias adjustment to stationarity
%%%%%%%%%%%% Note that the horizon h LP is a VARMA where the first h lags of y are zero (see e.g. Jorda and Kozicki (2011) equation 19). 
% This VARMA is inverted into it's corresponding VAR using the arma2ar function, and restricting bias to stationarity follows Kilian (1998)
% after setting up as a VAR(1). In order to avoid computational demands while avoiding practical issues of bias due to truncation,
% the VAR is truncated at lag T.

biaslp= biaslp/(T-arp-h); %estimate of bias
betahat_adj=betahat(2:end,:)-biaslp;

AR0b4=[zeros(q,q*(h-1)) betahat(2:end,:)'];

A=Beff(2:q+1,:,1:h-1,1);
C = permute(A,[1 3 2]);
C = reshape(C,[],size(A,2),1);

MA0=C';
AR0b4=mat2cell(AR0b4,q,q*ones(1,arp+h-1));
MA0=mat2cell(MA0,q,q*ones(1,h-1));
arb4 = arma2ar(AR0b4,MA0,T);
% arb4 = arma2ar(AR0b4,MA0);
arb4=cell2mat(arb4);
[qc, arpc]=size(arb4);
arpc=arpc/q;

Apbefore=arb4;
VAR_Companionbefore=[Apbefore; eye((arpc-1)*qc) zeros(((arpc-1)*qc),qc)];
 
 
if any(abs(eig(VAR_Companionbefore))>=1)
    Beff(2:end,:,h,1) = betahat(2:end,:) ;
else
    Beff(2:end,:,h,1) = betahat_adj ;
end
 

AR0=[zeros(q,q*(h-1)) Beff(2:end,:,h,1)'];

A=Beff(2:q+1,:,1:h-1,1);
C = permute(A,[1 3 2]);
C = reshape(C,[],size(A,2),1);

MA0=C';

AR0=mat2cell(AR0,q,q*ones(1,arp+h-1));
MA0=mat2cell(MA0,q,q*ones(1,h-1));
ar = arma2ar(AR0,MA0,T);
ar=cell2mat(ar);
[qc, arpc]=size(ar);
arpc=arpc/q;


Ap=ar;
VAR_Companion=[Ap; eye((arpc-1)*qc) zeros(((arpc-1)*qc),qc)];
 
john=0;
 
while all(abs(eig(VAR_Companionbefore))<1) && any(abs(eig(VAR_Companion))>=1)
john=john+1;
    Beff(2:end,:,h,1)=betahat(2:end,:)-(biaslp*(1-.01*john));
    
if john==100
        break
end  
    
AR0=[zeros(q,q*(h-1)) Beff(2:end,:,h,1)'];

A=Beff(2:q+1,:,1:h-1,1);
C = permute(A,[1 3 2]);
C = reshape(C,[],size(A,2),1);

MA0=C';

AR0=mat2cell(AR0,q,q*ones(1,arp+h-1));
MA0=mat2cell(MA0,q,q*ones(1,h-1));
ar = arma2ar(AR0,MA0,T);
ar=cell2mat(ar);
[qc, arpc]=size(ar);
arpc=arpc/q;

Ap=ar;
VAR_Companion=[Ap; eye((arpc-1)*qc) zeros(((arpc-1)*qc),qc)];
 
end
    end
%%%%%end of bias adjustment

  end
  
   
%%generate bootstrap samples for horizon 2 and greater
  for h=2:maxh
      
      for n=2:nstraps
          
        gammat2=zeros(p,q);

gammat1=(samecov*saveF_h1(:,1:end-h+1)*epshat(h:end,:,n))/(T-arp-h+1) ;  
h0=(samecov*saveF_h1(:,1:end-h+1)*epshat(1:end-h+1,:,n)/(T-arp-h+1));
    
if h>2
    for l=1:h-2
        gammat2=gammat2+(samecov*Qsome(:,:,h-l-1)'*h0*Beff(2:q+1,:,l,n));
    end
end
    
         l=h-1;
        gammat2=gammat2+(h0*Beff(2:q+1,:,l,n));
Beff(:,:,h,n)=Beff(:,:,h,1)+ gammat1+gammat2 ;

      end
      
  end
    

    if (full==0)
        Beff=Beff(2:q+1,:,:,:); %only keeping Wold IRFs
    
    elseif full==1
        %keeps all parameters not just Wold IRFs
    else 
        error('Error. Incorrect input for full. full must be 0 or 1')
    end

 
return
