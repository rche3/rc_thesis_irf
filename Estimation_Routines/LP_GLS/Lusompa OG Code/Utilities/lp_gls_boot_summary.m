function [CoverageLPtemp, LengthLPtemp, BiasLPtemp] = lp_gls_boot_summary(Beff,maxh,nstraps,confint,trueirf,q) 

%Calculates coverage, length, and bias for LP bootstrap estimator
%for the Monte Carlos (except for the LR restrictions case)

%Input:
    %Beff is the qxqx(maxh)xnstraps matrix from the lp_gls_boot function
    %maxh is the max horizon used in the lp_gls_boot function
    %nstraps is the number of bootstrap replications used in the lp_gls_boot function
    %confit in the size of the desired confidence interval (e.g. .95)
    %trueirf are the actual Wold irfs from the data generating process
    %q number of variables in the system used in the lp_gls_boot function

%Output:
    %CoverageLPtemp is a (q*q)x(maxh) matrix of coverage estimates (1 if confidence interval contains true estimate; 0 otherwise)
    %LengthLPtemp (q*q)x(maxh) matrix of length estimates
    %BiasLPtemp (q*q)x(maxh) matrix of bias estimates

CoverageLPtemp=zeros(maxh,q*q);
LengthLPtemp=zeros(maxh,q*q);
BiasLPtemp=zeros(maxh,q*q);

for h=1:maxh
  
    LPHolder=zeros(q*q,nstraps-1);
    
    for n=2:nstraps 

        if q>1
LPHolder(:,n-1)=reshape(Beff(:,:,h,n)',[],1);
        else
           LPHolder(:,n-1)=reshape(Beff(h,n)',[],1);
 
        end

    end

alpha=1-confint;
       lp_intervals= quantile(LPHolder',[(alpha/2) 1-(alpha/2)])';
                             
                      if q==1
                          
                                                LengthLPtemp(h,:)=lp_intervals(2)-lp_intervals(1);
          
          jack= (reshape(trueirf(h),[],1) >= lp_intervals(1)) .* (reshape(trueirf(h),[],1) <= lp_intervals(2));
          
                    CoverageLPtemp(h,:)=jack;
                    
                    Hold=mean(Beff(h,2:end),2);
                    
                                   BiasLPtemp(h,:)=Hold-trueirf(h); 
                   
                      else
                          
                                                LengthLPtemp(h,:)=(lp_intervals(:,2)-lp_intervals(:,1))';
                          
         jack= (reshape(trueirf(:,:,h),[],1) >= lp_intervals(:,1)) .* (reshape(trueirf(:,:,h),[],1) <= lp_intervals(:,2));

       
          CoverageLPtemp(h,:)=jack';
          
          Hold=mean(Beff(:,:,h,2:end),4);
          Hold=Hold';
          Hold=Hold-trueirf(:,:,h);
                         BiasLPtemp(h,:)=reshape(Hold,[],1)'; 
               
                       end
 
end

return