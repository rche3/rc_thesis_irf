function [CoverageLPHARtemp, LengthLPHARtemp, BiasLPHARtemp] = lp_har_summary_longrun(y,arp,maxh,maxhLR,confint,trueirf) 

%Calculates coverage, length, and bias for LP OLS estimator
%for the long run restriction Monte Carlo
    %Makes use of the symbolic toolbox

%Input:
    %y is the qxT matrix of variables
    %arp is the lag length
    %maxh is the maximum desired LP horizon
    %confint is desired confidence intervals for coverage calculation (e.g. .95)
    %trueirf are the true impulse responses functions

%Output:
    %CoverageLPHARtemp is a (q*q)x(maxh) matrix of coverage estimates (1 if confidence interval contains true estimate; 0 otherwise)
    %LengthLPHARtemp (q*q)x(maxh) matrix of length estimates
    %BiasLPHARtemp (q*q)x(maxh) matrix of bias estimates

Y=y-mean(y,2);  [q, T]=size(Y); 
saveY=Y; 
Teffective=T-arp;


alpha=1-confint;
upper=1-(alpha/2);

h=1;
    p=arp*q; 
F=zeros(q*arp,T-arp); 
for j=1:arp 
    F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-j)); 
end
 samecov=(F*F'/Teffective)\eye(p);

 yeff=saveY(:,arp+h:end)'; 


Teffective=size(yeff',2);
saveF=F;

    
[betahat,~,~,resids] = ols(yeff,F',0) ;

    Beff = zeros(p,q,maxh) ;
    Bigcov=zeros(p*q+(q*(q+1)/2),p*q+(q*(q+1)/2),maxh) ;
    Sigma_error=resids'*resids/Teffective;

    epshat=resids;

        Beff(:,:,h) = betahat ;


        
  dof = floor(0.41*T^(2/3));


Score=zeros(p*q+(q*(q+1)/2),Teffective);

  %calculating limiting variance for Sigma
  vecerror=zeros((q*(q+1)/2),Teffective);
vecsigma=reshape(Sigma_error,[],1);
            vecsigma(2)=[];
        for t=1:Teffective
            temp=reshape(epshat(t,:)'*epshat(t,:),[],1);
            temp(2)=[];

            vecerror(:,t)=temp-vecsigma;
        end


for t=1:Teffective
    Score(:,t)=[kron(samecov*saveF(:,t),eye(q))*resids(t,:)';vecerror(:,t)];
end


  % Cosine weights

  Lambda=zeros(p*q+(q*(q+1)/2),dof);
  for j=1:dof
       for t=1:Teffective
          Lambda(:,j)= Lambda(:,j)+ sqrt(2/Teffective)*Score(:,t)*cos(pi*j*(t-.5)/Teffective);
       end
       Bigcov(:,:,h)=Bigcov(:,:,h)+Lambda(:,j)*Lambda(:,j)'/dof;
  end



  

for h=2:maxhLR
  

    p=arp*q;
F=zeros(q*arp,T-arp-h+1); 
for j=1:arp 
    F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-h+1-j)); 
end


saveF=F;


yeff=saveY(:,arp+h:end)'; 


Teffective=size(yeff',2);

[betahat,~,~,resids] = ols(yeff,F',0) ;

    Beff(:,:,h) = betahat ;

Score=zeros(p*q+(q*(q+1)/2),Teffective);


for t=1:Teffective
    Score(:,t+h-1)=[kron(samecov*saveF(:,t),eye(q))*resids(t,:)';vecerror(:,t+h-1)];
end

if h<maxh+1 %only calculate for relevant horizons
  % Cosine weights

  Lambda=zeros(p*q+(q*(q+1)/2),dof);
  for j=1:dof
       for t=1:Teffective
          Lambda(:,j)= Lambda(:,j)+ sqrt(2/Teffective)*Score(:,t)*cos(pi*j*(t-.5)/Teffective);
       end
       Bigcov(:,:,h)=Bigcov(:,:,h)+Lambda(:,j)*Lambda(:,j)'/dof;
  end

end

end


%calculating coverage
CoverageLPHARtemp=zeros(maxh,q*q);
LengthLPHARtemp=zeros(maxh,q*q);
BiasLPHARtemp=zeros(maxh,q*q);


syms vechsigmaerror [(q*(q+1)/2) 1]


LRsum=eye(q);
       
for i=1:maxhLR
           LRsum=LRsum+Beff(1:q,:,i)';        
end

temp3=[vechsigmaerror1 vechsigmaerror2; vechsigmaerror2 vechsigmaerror3];
Hold=LRsum*temp3*LRsum';
P=LRsum\chol(Hold,'lower','nocheck');
vecP=reshape(P,[],1);

H=[];
for i=1:4
H=[H; transpose(gradient(vecP(i),vechsigmaerror))];
end

H=subs(H, {vechsigmaerror1,vechsigmaerror2,vechsigmaerror3}, {Sigma_error(1,1),Sigma_error(2,1),Sigma_error(2,2)});
H=double(H);  

P=subs(P, {vechsigmaerror1,vechsigmaerror2,vechsigmaerror3}, {Sigma_error(1,1),Sigma_error(2,1),Sigma_error(2,2)});
P=double(P);  


        
              

   Strucirf=zeros(q,q,maxh);
s_errors=zeros(q,q,maxh);
varcovirf=zeros(q*q,q*q,maxh);

for h=1:maxh
      Hold=Beff(1:q,:,h)';
Strucirf(:,:,h)=double(Hold*P);

    Hold=Beff(1:q,:,h)';
    

varcovirf(:,:,h)=(kron(P',eye(q))*Bigcov(1:q*q,1:q*q,h)*kron(P,eye(q))/Teffective)+(kron(eye(q),Hold)*H*Bigcov(end-(q*(q+1)/2)+1:end,end-(q*(q+1)/2)+1:end,h)*H'*kron(eye(q),Hold')/Teffective);

Hold=sqrt(diag(varcovirf(:,:,h)));
          s_errors(:,:,h)=reshape(Hold,[q,q]) ;

end

s_errors=squeeze(s_errors);


for h=1:maxh

            if q==1
                
                 lp_intervals=[Strucirf(h)-tinv(upper,dof)*s_errors(h)' Strucirf(h)+tinv(upper,dof)*s_errors(h)]; 
 
                                  LengthLPHARtemp(h,:)=lp_intervals(2)-lp_intervals(1);

                
                                    jack= (trueirf(h) >= lp_intervals(1)) .* (trueirf(h) <= lp_intervals(2)); 
          CoverageLPHARtemp(h,:)=jack';

               BiasLPHARtemp(h,:)=Strucirf(h)-trueirf(h); 
                
            else
                bottom=reshape(Strucirf(:,:,h),[],1)-(tinv(upper,dof)*reshape(s_errors(:,:,h),[],1));
                top=reshape(Strucirf(:,:,h),[],1)+(tinv(upper,dof)*reshape(s_errors(:,:,h),[],1));
                
                 lp_intervals=[bottom top];

 
                                  LengthLPHARtemp(h,:)=(lp_intervals(:,2)-lp_intervals(:,1))';

                    jack= (reshape(trueirf(:,:,h),[],1) >= lp_intervals(:,1)) .* (reshape(trueirf(:,:,h),[],1) <= lp_intervals(:,2)); 
          CoverageLPHARtemp(h,:)=jack';

               BiasLPHARtemp(h,:)=reshape(Strucirf(:,:,h),[],1)-reshape(trueirf(:,:,h),[],1); 
            
            end

end


return