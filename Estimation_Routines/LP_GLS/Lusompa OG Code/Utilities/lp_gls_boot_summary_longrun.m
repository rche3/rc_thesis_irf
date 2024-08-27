function [CoverageLPtemp, LengthLPtemp, BiasLPtemp] = lp_gls_boot_summary_longrun(Beff,Sigma_hold,maxh,maxhLR,nstraps,confint,trueirf,q) 

%Calculates coverage, length, and bias for LP bootstrap estimator
%for the long run restrictions Monte Carlo

%Input:
    %Beff are LP estimates output from the lp_gls_boot function
    %Sigma_hold the residual covariance matrix output from the lp_gls_analytical_longrun function
    %maxh is maximum LP horizon used in the lp_gls_boot function
    %nstraps is the number of bootstrap replications used in the lp_gls_boot function
    %confint is desired confidence intervals for coverage calculation (e.g. .95)
    %trueirf are the true impulse responses functions
    %q is the number of variables used in the lp_gls_boot function

%Output:
    %CoverageLPtemp is a (q*q)x(maxh) matrix of coverage estimates (1 if confidence interval contains true estimate; 0 otherwise)
    %LengthLPtemp (q*q)x(maxh) matrix of length estimates
    %BiasLPtemp (q*q)x(maxh) matrix of bias estimates

CoverageLPtemp=zeros(maxh,q*q);
LengthLPtemp=zeros(maxh,q*q);
BiasLPtemp=zeros(maxh,q*q);


%imposing long run restrictions 
Cholder=zeros(q,q,nstraps-1);


LRsum=repmat(eye(q),1,nstraps-1);
LRsum=reshape(LRsum,[q q (nstraps-1)]);
       
for n=2:nstraps 
for i=1:maxhLR
           LRsum(:,:,n-1)=LRsum(:,:,n-1)+ Beff(:,:,i,n)';
end

Hold=LRsum(:,:,n-1)*Sigma_hold(:,:,n)*LRsum(:,:,n-1)';
Cholder(:,:,n-1)=LRsum(:,:,n-1)\chol(Hold,'lower');
end
%end of imposing long run restrictions 

%%calculates coverage, length, and bias 
for h=1:maxh

    LPHolder=zeros(q*q,nstraps-1);
    
    for n=2:nstraps 

LPHolder(:,n-1)=reshape(Beff(:,:,h,n)'*Cholder(:,:,n-1),[],1);

    end

alpha=1-confint;
       lp_intervals= quantile(LPHolder',[(alpha/2) 1-(alpha/2)])';
       
                         
                                                LengthLPtemp(h,:)=(lp_intervals(:,2)-lp_intervals(:,1))';
                      
         jack= (reshape(trueirf(:,:,h),[],1) >= lp_intervals(:,1)) .* (reshape(trueirf(:,:,h),[],1) <= lp_intervals(:,2));
    
          CoverageLPtemp(h,:)=jack';
          
          Hold=mean(LPHolder,2);
          Hold=reshape(Hold,q,q);
          Hold=Hold-trueirf(:,:,h);
                         BiasLPtemp(h,:)=reshape(Hold,[],1)'; 



end

return