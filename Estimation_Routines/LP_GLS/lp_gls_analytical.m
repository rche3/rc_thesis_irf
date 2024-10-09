function [Beff, s_errors] = lp_gls_analytical(y,arp,maxh,badj) 

%Analytical LP GLS Code from Lusompa (2023)

%Input:
    %y is the qxT matrix of variables
    %arp is the lag length
    %maxh is the maximum desired LP horizon
    %badj equals 1 for bias adjustment and 0 otherwise

%Output:
    %Beff is the qxqxmaxh matrix of the wold impulse responses bootstrap replications. For
        %example Beff(:,:,h)' gives the horizon h estimate of the wold impulse
    %s_errors qxqxmaxh matrix of standard errors corresponding to the Beff matrix

Y=y; [q, T]=size(Y); 
  

%population estimates for h=1
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

epshat=resids;

    Beff = zeros(p,q,maxh) ;
    s_errors=zeros(p,q,maxh) ;


    Beff(:,:,h) = betahat ;
        Sigma_hold=resids'*resids/Teffective;
        %bias adjusment
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

    
%%%%%%%%%%%%
biaslp= biaslp./Teffective;  %estimate of bias
betahat_adj=betahat(2:end,:)-biaslp;
 
 
Apbefore=betahat(2:end,:)';
VAR_Companionbefore=[Apbefore; eye((arp-1)*q) zeros(((arp-1)*q),q)];
 
%restricting bias to stationarity 
if any(abs(eig(VAR_Companionbefore))>=1)
    Beff(2:end,:,h) = betahat(2:end,:) ;
else
    Beff(2:end,:,h) = betahat_adj ;
end
 
Ap=Beff(2:end,:,h)';
VAR_Companion=[Ap; eye((arp-1)*q) zeros(((arp-1)*q),q)];
 
john=0;
 
while all(abs(eig(VAR_Companionbefore))<1) && any(abs(eig(VAR_Companion))>=1)
john=john+1;
    Beff(2:end,:,h)=betahat(2:end,:)-(biaslp*(1-.01*john));
    
if john==100
        break
end  
    
    Ap=Beff(2:end,:,h)';
VAR_Companion=[Ap; eye((arp-1)*q) zeros(((arp-1)*q),q)];

end
%%%%% end of bias adjustment
 
      epshat = yeff-F'*Beff(:,:,1);  %population residuals for horizon 0
      
      epshat=epshat-mean(epshat);
end

 
%Calculating SE for horizon 1    
    saveF_h1=F;
    samecov=(saveF_h1*saveF_h1'/Teffective)\eye(p);
    
    Qsome=zeros(p,p,maxh-1);
 
    for l=1:maxh-1
        Qsome(:,:,l)=saveF_h1(:,1+l:end)*(saveF_h1(:,1:end-l)'/(Teffective)); 
    end
       
  Hold=zeros(p*q,p*q);
  
  for t=1:T-arp
       Hold=Hold+(kron(eye(q),samecov*saveF_h1(:,t))*(epshat(t,:)'*epshat(t,:))*kron(eye(q),samecov*saveF_h1(:,t))');
  end
    
    Hold=Hold/((T-arp)*(T-arp));
    Hold=sqrt(diag(Hold));
          s_errors(:,:,h)=reshape(Hold,[p,q]) ;
           
end
%end of Calculating SE for horizon 1    
%end of population estimates for horizon 1

 
 
 
%population estimate for h>=2 
  for h = 2:maxh 
      
          epscorrect = zeros(T-arp-h+1,q) ;  
          
    %gls correction for population estimate
    for hsub = 1:h-1 
      epscorrect = epscorrect +  epshat(h-hsub:T-arp-hsub,:)*Beff(2:q+1,:,hsub) ;
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
F=[ones(1,Teffective);saveF]; % Design Matrix
 
 
    %population estimate
    [betahat,~,~,resids] = ols(yeff,F',0) ;
    
    Beff(:,:,h) = betahat ;
 

%bias adjusment for horizon h
    if badj==1

        [A,SIGMA,~,~,~]=olsvarc(y(:,1:end-h)',arp);

phitwid=A;
omegaUtwid=SIGMA;
nZtwid=size(SIGMA,1);
nk=nZtwid;
PX=eye(nZtwid);

        biaslp=zeros(p-1,q);


        %calculating bias for each equation
for dim=1:q

EetaZtwid0=(F(2:end,arp+1:end)*resids(1:end-arp,dim))'/(Teffective-arp);

[vbias] = proc_vb_ma0(nZtwid, phitwid, omegaUtwid, nk, PX, EetaZtwid0);  %bias adjustment function

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

A=Beff(2:q+1,:,1:h-1);
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
    Beff(2:end,:,h) = betahat(2:end,:) ;
else
    Beff(2:end,:,h) = betahat_adj ;
end


AR0=[zeros(q,q*(h-1)) Beff(2:end,:,h)'];

A=Beff(2:q+1,:,1:h-1);
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
    Beff(2:end,:,h)=betahat(2:end,:)-(biaslp*(1-.01*john));
    
if john==100
        break
end  
    
AR0=[zeros(q,q*(h-1)) Beff(2:end,:,h)'];

A=Beff(2:q+1,:,1:h-1);
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


%Calculating SE for horizon h
  Hold=zeros(p*q,p*q);  
  sk=0;
  
  if h>2
    for l=1:h-2
        sk=sk+kron((Beff(2:q+1,:,l)'),samecov*Qsome(:,:,h-l-1)');
    end
  end    
         l=h-1;
        sk=sk+kron((Beff(2:q+1,:,l)'),eye(p));
               
  for t=1:T-arp-h+1


      Hold=Hold+(kron(eye(q),samecov*saveF_h1(:,t))*(epshat(t+h-1,:)'*epshat(t+h-1,:))*kron(eye(q),samecov*saveF_h1(:,t))');
      
            Hold=Hold+(sk*kron(eye(q),samecov*saveF_h1(:,t))*(epshat(t,:)'*epshat(t,:))*kron(eye(q),samecov*saveF_h1(:,t))'*sk');
            
                  Hold=Hold+(kron(eye(q),samecov*saveF_h1(:,t))*(epshat(t+h-1,:)'*epshat(t+h-1,:))*kron(eye(q),samecov*saveF_h1(:,t+h-1))'*sk');
 
                        Hold=Hold+(sk*kron(eye(q),samecov*saveF_h1(:,t+h-1))*(epshat(t+h-1,:)'*epshat(t+h-1,:))*kron(eye(q),samecov*saveF_h1(:,t))');

  end
     
    Hold=Hold/((T-arp-h+1)*(T-arp-h+1));
    Hold=sqrt(diag(Hold));
          s_errors(:,:,h)=reshape(Hold,[p,q]) ;
%end of calculating standard errors for horizon h
  end
%end of population estimate for h>=2


  %only keeping Wold IRFs and corresponding SE
  Beff=Beff(2:q+1,:,:); 
    s_errors=s_errors(2:q+1,:,:);

return
