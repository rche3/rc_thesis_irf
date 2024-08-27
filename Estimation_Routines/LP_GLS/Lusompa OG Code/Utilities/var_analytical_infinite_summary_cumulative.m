function [CoverageVAR_analytictemp, LengthVAR_analytictemp, BiasVAR_analytictemp] = var_analytical_infinite_summary_cumulative(y,arp,maxh,trueirf,confint) 

%Calculates coverage, length, and bias in Monte Carlos for infinite order VAR following Lutkepohl (1990) Prop 1
%and Lutkepohl 2005 prop 15.4. Calculate summary statistics in long run restriction case.
%VAR estimates and estimate covariance matrices are calculated before coverage

%Input:
    %y is the qxT matrix of variables
    %arp is the lag length
    %maxh is maximum horizon to be calculated for coverage, length, and bias
    %trueirf are the true impulse responses functions
    %confint is desired confidence intervals for coverage calculation (e.g. .95)

%Output:
    %CoverageVARtemp is a (q*q)x(maxh) matrix of coverage estimates (1 if confidence interval contains true estimate; 0 otherwise)
    %LengthVARtemp (q*q)x(maxh) matrix of length estimates
    %BiasVARtemp (q*q)x(maxh) matrix of bias estimates
    
    Y=y-mean(y,2); [q, T]=size(Y); 

             p=arp*q; 
             F=zeros(q*arp,T-arp); 
             for j=1:arp
                F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-j)); 
             end
             Y(:,1:arp)=[]; Teffective=size(Y,2);              
             

% estimate rho:
  
  [betah,~,~,epshat]= ols(Y',F',0) ;
      
    samecov=(F*F'/Teffective)\eye(p); %Gamma inverse in the paper
    
Sigma_error=epshat'*epshat/(Teffective);
Sigma_error_inv=Sigma_error\eye(q);
Sigma_alpha=kron(samecov,Sigma_error);
    

        Ap=betah';
VAR_Companion=[Ap; eye((arp-1)*q) zeros(((arp-1)*q),q)];

  
  Beff=zeros(q,q,maxh);
    Beff_alt=zeros(q,q,maxh+1);
    Beff_alt(:,:,1)=eye(q);

varcovirf=zeros(q*q,q*q,maxh);
s_errors=zeros(q,q,maxh);

for h=1:maxh
      Hold=VAR_Companion^h;

      if h==1
Beff(:,:,h)=Hold(1:q,1:q);
      else
Beff(:,:,h)=Beff(:,:,h-1)+Hold(1:q,1:q);
      end

      Beff_alt(:,:,h+1)=Hold(1:q,1:q);

G=0;

for n=2:h+1
    for m=2:h+1
        for i=1:m
            if n-m+i<1
            G=G+Beff_alt(:,:,i)*Sigma_error*zeros(q);
            else
            G=G+Beff_alt(:,:,i)*Sigma_error*Beff_alt(:,:,n-m+i)';
            end
        end
    end
end


varcovirf(:,:,h)=kron(Sigma_error_inv,G)/Teffective;

Hold=sqrt(diag(varcovirf(:,:,h)));
          s_errors(:,:,h)=reshape(Hold,[q,q]) ;

end

Beff=squeeze(Beff);
s_errors=squeeze(s_errors);





%%%%% Calculating coverage, length, and bias
  CoverageVAR_analytictemp=zeros(maxh,q*q);
LengthVAR_analytictemp=zeros(maxh,q*q);
BiasVAR_analytictemp=zeros(maxh,q*q);




for h=1:maxh
    
alpha=1-confint;
upper=1-(alpha/2);


            if q==1
                
                 var_intervals=[Beff(h)-tinv(upper,T-(arp*q)-1)*s_errors(h) Beff(h)+tinv(upper,T-(arp*q)-1)*s_errors(h)]; 
 
                                  LengthVAR_analytictemp(h,:)=var_intervals(2)-var_intervals(1);

                
                                    jack= (trueirf(h) >= var_intervals(1)) .* (trueirf(h) <= var_intervals(2)); 
          CoverageVAR_analytictemp(h,:)=jack';

               BiasVAR_analytictemp(h,:)=Beff(h)-trueirf(h); 
                
            else
                bottom=reshape(Beff(:,:,h),[],1)-(tinv(upper,T-(arp*q)-1)*reshape(s_errors(:,:,h),[],1));
                top=reshape(Beff(:,:,h),[],1)+(tinv(upper,T-(arp*q)-1)*reshape(s_errors(:,:,h),[],1));
                
                 var_intervals=[bottom top];

 
                                  LengthVAR_analytictemp(h,:)=(var_intervals(:,2)-var_intervals(:,1))';

                    jack= (reshape(trueirf(:,:,h),[],1) >= var_intervals(:,1)) .* (reshape(trueirf(:,:,h),[],1) <= var_intervals(:,2)); 
          CoverageVAR_analytictemp(h,:)=jack';

               BiasVAR_analytictemp(h,:)=reshape(Beff(:,:,h),[],1)-reshape(trueirf(:,:,h),[],1); 
            
            end
            
end


 


return