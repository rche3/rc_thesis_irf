function [CoverageLPtemp, LengthLPtemp, BiasLPtemp] = lp_gls_analytical_summary_longrun(Beff,Bigcov,Sigma_hold,epshat,maxh,maxhLR,confint,trueirf,T,arp,q)

%Calculates coverage, length, and bias for LP analytical estimator
%for the long run restrictions Monte Carlo
%Makes use of the symbolic toolbox

%Input:
    %Beff are LP estimates output from the lp_gls_analytical_longrun function
    %Bigcov is LP covariance matrix output from the lp_gls_analytical_longrun function
    %Sigma_hold the residual covariance matrix output from the lp_gls_analytical_longrun function
    %epshat is the matrix of residuals output from the lp_gls_analytical_longrun function
    %maxh is the maximum desired LP horizon (should be same one used for the lp_gls_analytical_longrun function)
    %confint is desired confidence intervals for coverage calculation (e.g. .95)
    %trueirf are the true impulse responses functions
    %T is the sample size of data used in lp_gls_analytical_longrun function
    %arp is the lag length used in lp_gls_analytical_longrun function
    %q is the number of variables used in the lp_gls_analytical_longrun function

%Output:
    %CoverageLPtemp is a (q*q)x(maxh) matrix of coverage estimates (1 if confidence interval contains true estimate; 0 otherwise)
    %LengthLPtemp (q*q)x(maxh) matrix of length estimates
    %BiasLPtemp (q*q)x(maxh) matrix of bias estimates


CoverageLPtemp=zeros(maxh,q*q);
LengthLPtemp=zeros(maxh,q*q);
BiasLPtemp=zeros(maxh,q*q);

Teffective=T-arp;

%imposing long run restrictions and calculating SE
%Is simply the LP application of Lutkepohl (1988) "Asymptotic Distribution of the Moving
%Average Coefficients of an Estimated Vector Autoregressive Process".

syms vechsigmaerror [3 1]

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

H=subs(H, {vechsigmaerror1,vechsigmaerror2,vechsigmaerror3}, {Sigma_hold(1,1),Sigma_hold(2,1),Sigma_hold(2,2)});
H=double(H);  

P=subs(P, {vechsigmaerror1,vechsigmaerror2,vechsigmaerror3}, {Sigma_hold(1,1),Sigma_hold(2,1),Sigma_hold(2,2)});
P=double(P);  

epshatsq=zeros(3,3,Teffective);
        for t=1:Teffective
            vecerror=reshape(epshat(t,:)'*epshat(t,:),[],1);
            vecerror(2)=[];

            temp=vecerror*vecerror';
       epshatsq(:,:,t)=temp;
        end
        
            vecsigma=reshape(Sigma_hold,[],1);
            vecsigma(2)=[];
 Sigma_sigma=mean(epshatsq,3)-vecsigma*vecsigma';

 p=arp*q;
  

   Strucirf=zeros(q,q,maxh);
s_errors=zeros(q,q,maxh);
lpcovirf=zeros(q*q,q*q,maxh);


for h=1:maxh
      Hold=Beff(1:q,:,h)';
Strucirf(:,:,h)=double(Hold*P);

    
lpcovirf(:,:,h)=(kron(P',eye(q))*Bigcov(1:q*q,1:q*q,h)*kron(P,eye(q))/Teffective)+(kron(eye(q),Hold)*H*Sigma_sigma*H'*kron(eye(q),Hold')/Teffective);

Hold=sqrt(diag(lpcovirf(:,:,h)));
          s_errors(:,:,h)=reshape(Hold,[q,q]) ;

end
%end of imposing long run restrictions and calculating SE

s_errors=squeeze(s_errors);


%%calculates coverage, length, and bias 
for h=1:maxh
    
alpha=1-confint;
upper=1-(alpha/2);

            if q==1
                
                 lp_intervals=[Strucirf(h)-tinv(upper,T-(arp*q)-h)*s_errors(h) Strucirf(h)+tinv(upper,T-(arp*q)-h)*s_errors(h)]; 
 
                                  LengthLPtemp(h,:)=lp_intervals(2)-lp_intervals(1);

                
                                    jack= (trueirf(h) >= lp_intervals(1)) .* (trueirf(h) <= lp_intervals(2)); 
          CoverageLPtemp(h,:)=jack';

               BiasLPtemp(h,:)=Strucirf(h)-trueirf(h); 
                
            else
                bottom=reshape(Strucirf(:,:,h),[],1)-(tinv(upper,T-(arp*q)-h)*reshape(s_errors(:,:,h),[],1));
                top=reshape(Strucirf(:,:,h),[],1)+(tinv(upper,T-(arp*q)-h)*reshape(s_errors(:,:,h),[],1));
                
                 lp_intervals=[bottom top];

 
                                  LengthLPtemp(h,:)=(lp_intervals(:,2)-lp_intervals(:,1))';

                    jack= (reshape(trueirf(:,:,h),[],1) >= lp_intervals(:,1)) .* (reshape(trueirf(:,:,h),[],1) <= lp_intervals(:,2)); 
          CoverageLPtemp(h,:)=jack';

               BiasLPtemp(h,:)=reshape(Strucirf(:,:,h),[],1)-reshape(trueirf(:,:,h),[],1); 
            
            end



end

return