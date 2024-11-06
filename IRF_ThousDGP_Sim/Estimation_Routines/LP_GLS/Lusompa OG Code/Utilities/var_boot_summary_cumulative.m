function [CoverageVARtemp, LengthVARtemp, BiasVARtemp] = var_boot_summary_cumulative(rhohat,Sigma_hold,maxh,nstraps,confint,arp,q,trueirf) 

%Calculates coverage, length, and bias for LP bootstrap estimator
%for the long run restrictions Monte Carlo

%Input:
    %rhohat are VAR estimates output from the var_boot function
    %Sigma_hold the residual covariance matrix output from the var_boot function
    %maxh is maximum horizon to be used when calculating coverage, length, and bias
    %nstraps is the number of bootstrap replications used in the var_boot function
    %confint is desired confidence intervals for coverage calculation (e.g. .95)
    %arp is the lag length used in var_boot function
    %q is the number of variables used in the var_boot function
    %trueirf are the true impulse responses functions

%Output:
    %CoverageVARtemp is a (q*q)x(maxh) matrix of coverage estimates (1 if confidence interval contains true estimate; 0 otherwise)
    %LengthVARtemp (q*q)x(maxh) matrix of length estimates
    %BiasVARtemp (q*q)x(maxh) matrix of bias estimates

CoverageVARtemp=zeros(maxh,q*q);
LengthVARtemp=zeros(maxh,q*q);
BiasVARtemp=zeros(maxh,q*q);


for h=1:maxh
  
    VARHolder=zeros(q*q,nstraps-1);
    
    for n=2:nstraps 
    Ap=rhohat(2:end,:,n)';
VAR_Companion=[Ap; eye((arp-1)*q) zeros(((arp-1)*q),q)];

Holder=zeros(size(VAR_Companion,1));
for i=1:h
Holder=Holder+(VAR_Companion^i);
end

VARHolder(:,n-1)=reshape(Holder(1:q,1:q),[],1);


    end

alpha=1-confint;

       var_intervals= quantile(VARHolder',[(alpha/2) 1-(alpha/2)])';
       
                                 
                               if q==1
                                   
                                   LengthVARtemp(h,:)=(var_intervals(2)-var_intervals(1))';
 
          jack= (trueirf(h) >= var_intervals(1)) .* (trueirf(h) <= var_intervals(2));
          CoverageVARtemp(h,:)=jack';
          
          
    Ap=rhohat(2:end,:,1)';
VAR_Companion=[Ap; eye((arp-1)*q) zeros(((arp-1)*q),q)];

Holder=VAR_Companion^h;
Holder=Holder(1:q,1:q);
                        
          BiasVARtemp(h,:)=Holder-trueirf(h); 
          
                               else
                                   LengthVARtemp(h,:)=(var_intervals(:,2)-var_intervals(:,1))';
                                   
                                             jack= (reshape(trueirf(:,:,h),[],1) >= var_intervals(:,1)) .* (reshape(trueirf(:,:,h),[],1) <= var_intervals(:,2));
          CoverageVARtemp(h,:)=jack';
          
          
    Ap=rhohat(2:end,:,1)';
VAR_Companion=[Ap; eye((arp-1)*q) zeros(((arp-1)*q),q)];

Holder=zeros(size(VAR_Companion,1));
for i=1:h
Holder=Holder+(VAR_Companion^i);
end                        
          BiasVARtemp(h,:)=reshape(Holder(1:q,1:q),[],1)-reshape(trueirf(:,:,h),[],1); 
                                   
                               end

end

return