function [CoverageLPtemp, LengthLPtemp, BiasLPtemp] = lp_gls_analytical_summary(Beff,s_errors,maxh,confint,trueirf,T,arp,q) 

%Calculates coverage, length, and bias for LP analytical estimator 
%for the Monte Carlos (except for the LR restrictions case)

%Input:
    %Beff are LP estimates output from the lp_gls_analytical function
    %s_errors are the standard errors output from the lp_gls_analytical function
    %maxh is maximum LP horizon used in the lp_gls_analytical function
    %confint is desired confidence intervals for coverage calculation (e.g. .95)
    %trueirf are the true impulse responses functions
    %T is the sample size of data used in lp_gls_analytical function
    %arp is the lag length used in lp_gls_analytical function
    %q is the number of variables used in the lp_gls_analytical function

%Output:
    %CoverageLPtemp is a (q*q)x(maxh) matrix of coverage estimates (1 if confidence interval contains true estimate; 0 otherwise)
    %LengthLPtemp (q*q)x(maxh) matrix of length estimates
    %BiasLPtemp (q*q)x(maxh) matrix of bias estimates


CoverageLPtemp=zeros(maxh,q*q);
LengthLPtemp=zeros(maxh,q*q);
BiasLPtemp=zeros(maxh,q*q);

%%calculates coverage, length, and bias 

for h=1:maxh
    
alpha=1-confint;
upper=1-(alpha/2);

            if q==1
                
                 lp_intervals=[Beff(h)-tinv(upper,T-(arp*q)-h)*s_errors(h)' Beff(h)+tinv(upper,T-(arp*q)-h)*s_errors(h)]; 
 
                                  LengthLPtemp(h,:)=lp_intervals(2)-lp_intervals(1);
             
                                    jack= (trueirf(h) >= lp_intervals(1)) .* (trueirf(h) <= lp_intervals(2)); 
          CoverageLPtemp(h,:)=jack';

               BiasLPtemp(h,:)=Beff(h)-trueirf(h); 
                
            else
                bottom=reshape(Beff(:,:,h)',[],1)-(tinv(upper,T-(arp*q)-h)*reshape(s_errors(:,:,h)',[],1));
                top=reshape(Beff(:,:,h)',[],1)+(tinv(upper,T-(arp*q)-h)*reshape(s_errors(:,:,h)',[],1));
                
                 lp_intervals=[bottom top];

                                  LengthLPtemp(h,:)=(lp_intervals(:,2)-lp_intervals(:,1))';

                    jack= (reshape(trueirf(:,:,h),[],1) >= lp_intervals(:,1)) .* (reshape(trueirf(:,:,h),[],1) <= lp_intervals(:,2)); 
          CoverageLPtemp(h,:)=jack';

               BiasLPtemp(h,:)=reshape(Beff(:,:,h)',[],1)-reshape(trueirf(:,:,h),[],1); 
            
            end


end

return