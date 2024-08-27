function [Beff, Bigcov,Sigma_hold,epshat] = lp_gls_analytical_longrun(y,arp,maxh,badj) 

%Calculates LP estimates to be used in the long run restrictions Monte
%Carlo

%Input:
    %y is the qxT matrix of variables
    %arp is the lag length
    %maxh is the maximum desired LP horizon
    %badj equals 1 for bias adjustment and 0 otherwise

%Output:
    %Beff is the pxqx(maxh) matrix of the LP estimates for horizon h.
    %Bigcov is a (p*q)x(p*q) by maxh covariance matrix for Beff' (note the transpose after Beff)
    %Sigma_hold qxq covariance matrix of the VAR (horizon 1 LP) residuals
    %epshat qxT matrix of the VAR (horizon 1 LP) residuals




%demeaning the data so intercept is not estimated for any horizon
Y=y-mean(y,2); [q, T]=size(Y); 
  
 
%population estimates for h=1 
for h = 1 
      
saveY=Y; 
   
    p=arp*q; 
F=zeros(q*arp,T-arp-h+1); 
for j=1:arp 
    F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-h+1-j)); 
end
 
saveF=F; % Design Matrix
 
 
yeff=saveY(:,arp+h:end)'; %Dependent variable
  
Teffective=size(yeff',2);
  
[betahat,~,~,resids] = ols(yeff,F',0) ;

epshat=resids;
 
    Beff = zeros(p,q,maxh) ;
    Bigcov=zeros(p*q,p*q,maxh) ;

 
    Beff(:,:,h) = betahat ;
        Sigma_hold(:,:)=resids'*resids/Teffective;


                %bias adjusment
if badj==1

    biaslp=zeros(p,q);

[A,SIGMA,~,~,~]=olsvarc(Y(:,1:end-h)',arp);

phitwid=A;
omegaUtwid=SIGMA;
nZtwid=size(SIGMA,1);
nk=nZtwid;
PX=eye(nZtwid);

%calculating bias for each equation
for dim=1:q

EetaZtwid0=(F(:,arp+1:end)*resids(1:end-arp,dim))'/(Teffective-arp);

[vbias] = proc_vb_ma0(nZtwid, phitwid, omegaUtwid, nk, PX, EetaZtwid0); %bias adjustment function

biaslp(:,dim)=vbias;

end
    
biaslp= biaslp./Teffective; %estimate of bias
betahat_adj=betahat-biaslp;
 
 
Apbefore=betahat';
VAR_Companionbefore=[Apbefore; eye((arp-1)*q) zeros(((arp-1)*q),q)];
 
%restricting bias to stationarity  
if any(abs(eig(VAR_Companionbefore))>=1)
    Beff(:,:,h) = betahat ;
else
    Beff(:,:,h) = betahat_adj ;
end
 
Ap=Beff(:,:,h)';
VAR_Companion=[Ap; eye((arp-1)*q) zeros(((arp-1)*q),q)];
 
john=0;
 
while all(abs(eig(VAR_Companionbefore))<1) && any(abs(eig(VAR_Companion))>=1)
john=john+1;
    Beff(:,:,h)=betahat-(biaslp*(1-.01*john));
    
if john==100
        break
end  
    
    Ap=Beff(:,:,h)';
VAR_Companion=[Ap; eye((arp-1)*q) zeros(((arp-1)*q),q)];
 
end
%%%%%end of bias adjustment
 
 
      epshat = yeff-F'*Beff(:,:,1);  %population residuals for horizon 1
      
      epshat=epshat-mean(epshat);
end

 
%Calculating covariance matrix of LP estimates for horizon 1    

    saveF_h1=F;
    samecov=(saveF_h1*saveF_h1'/Teffective)\eye(p);
    
    Qsome=zeros(p,p,maxh-1);
 
    for l=1:maxh-1
        Qsome(:,:,l)=saveF_h1(:,1+l:end)*(saveF_h1(:,1:end-l)'/(Teffective)); %could also divide by (T-arp) for "unbiased" degrees of freedom correction
    end
    
  Hold=zeros(p*q,p*q);
  
  for t=1:T-arp

                  Hold=Hold+(kron(eye(q),samecov*saveF_h1(:,t))*(epshat(t,:)'*epshat(t,:))*kron(eye(q),samecov*saveF_h1(:,t))');

  end
    
    Hold=Hold/(T-arp);
          Bigcov(:,:,h)=Hold ;
          
end
%end of Calculating SE for horizon 1    
%end of population estimates for horizon 1

 
 %population estimate for h>=2
  for h = 2:maxh 
      
          epscorrect = zeros(T-arp-h+1,q) ;  %Since you write h0 as h1 here you need to add 1 hence T-arp-h+1 instead of T-arp-h
          
    %gls correction for population estimate
    for hsub = 1:h-1 
      epscorrect = epscorrect +  epshat(h-hsub:T-arp-hsub,:)*Beff(1:q,:,hsub) ;
    end 
    
     
yeff=saveY(:,arp+h:end)'-epscorrect; %Dependent Variable

    p=arp*q; 
F=zeros(q*arp,T-arp-h+1); 
for j=1:arp 
    F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-h+1-j));  
end % Design Matrix
  
Teffective=size(yeff',2);
  
    %population estimate
    [betahat,~,~,resids] = ols(yeff,F',0) ;
    
    Beff(:,:,h) = betahat ;
 
        

%bias adjusment for horizon h
    if badj==1

        [A,SIGMA,~,~,~]=olsvarc(Y(:,1:end-h)',arp);

phitwid=A;
omegaUtwid=SIGMA;
nZtwid=size(SIGMA,1);
nk=nZtwid;
PX=eye(nZtwid);

biaslp=zeros(p-1,q);

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

biaslp= biaslp/(T-arp-h);  %estimate of bias
betahat_adj=betahat-biaslp;

AR0b4=[zeros(q,q*(h-1)) betahat'];

A=Beff(1:q,:,1:h-1);
C = permute(A,[1 3 2]);
C = reshape(C,[],size(A,2),1);

MA0=C';
AR0b4=mat2cell(AR0b4,q,q*ones(1,arp+h-1));
MA0=mat2cell(MA0,q,q*ones(1,h-1));
arb4 = arma2ar(AR0b4,MA0,T);
arb4=cell2mat(arb4);
[qc, arpc]=size(arb4);
arpc=arpc/q;

Apbefore=arb4;
VAR_Companionbefore=[Apbefore; eye((arpc-1)*qc) zeros(((arpc-1)*qc),qc)];
 

if any(abs(eig(VAR_Companionbefore))>=1)
    Beff(:,:,h) = betahat ;
else
    Beff(:,:,h) = betahat_adj ;
end
 
AR0=[zeros(q,q*(h-1)) Beff(2:end,:,h)'];

A=Beff(1:q,:,1:h-1);
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
    Beff(:,:,h)=betahat(:,:)-(biaslp*(1-.01*john));
    
if john==100
        break
end  
    
AR0=[zeros(q,q*(h-1)) Beff(:,:,h)'];

A=Beff(1:q,:,1:h-1);
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
%%%%%end of bias adjustment for horizon h


%Calculating covariance matrix for horizon h LP
  Hold=zeros(p*q,p*q);
  
  sk=0;
  
  if h>2
    for l=1:h-2
        sk=sk+kron(samecov*Qsome(:,:,h-l-1)',(Beff(1:q,:,l)'));
    end
  end

         l=h-1;
        sk=sk+kron(eye(p),(Beff(1:q,:,l)'));
          
  for t=1:T-arp-h+1

      Hold=Hold+(kron(samecov*saveF_h1(:,t),eye(q))*(epshat(t+h-1,:)'*epshat(t+h-1,:))*kron(samecov*saveF_h1(:,t),eye(q))');
      
            Hold=Hold+(sk*kron(samecov*saveF_h1(:,t),eye(q))*(epshat(t,:)'*epshat(t,:))*kron(samecov*saveF_h1(:,t),eye(q))'*sk');
            
                  Hold=Hold+(kron(samecov*saveF_h1(:,t),eye(q))*(epshat(t+h-1,:)'*epshat(t+h-1,:))*kron(samecov*saveF_h1(:,t+h-1),eye(q))'*sk');
             
                        Hold=Hold+(sk*kron(samecov*saveF_h1(:,t+h-1),eye(q))*(epshat(t+h-1,:)'*epshat(t+h-1,:))*kron(samecov*saveF_h1(:,t),eye(q))');

  end
  
    
    Hold=Hold/(T-arp-h+1);
          Bigcov(:,:,h)=Hold ;

  end
  %end of population estimate for h>=2
 

return
