function [CoverageVAR_analytictemp, LengthVAR_analytictemp, BiasVAR_analytictemp] = var_analytical_infinite_summary_longrun(y,arp,maxh,trueirf,confint) 

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
      

    
Sigma_error=epshat'*epshat/(Teffective);
% Sigma_error_inv=Sigma_error\eye(q);
    

        Ap=betah';
VAR_Companion=[Ap; eye((arp-1)*q) zeros(((arp-1)*q),q)];

syms vechsigmaerror [3 1]


LRsum=eye(q);
       
for i=1:arp
           LRsum=LRsum- Ap(:,(i-1)*q+(1:q),:);        
end

temp=LRsum\eye(q);
temp3=[vechsigmaerror1 vechsigmaerror2; vechsigmaerror2 vechsigmaerror3];
Hold=temp*temp3*temp';
P=LRsum*chol(Hold,'lower','nocheck');
vecP=reshape(P,[],1);

H=[];
for i=1:4
H=[H; transpose(gradient(vecP(i),vechsigmaerror))];
end

H=subs(H, {vechsigmaerror1,vechsigmaerror2,vechsigmaerror3}, {Sigma_error(1,1),Sigma_error(2,1),Sigma_error(2,2)});
H=double(H);  

P=subs(P, {vechsigmaerror1,vechsigmaerror2,vechsigmaerror3}, {Sigma_error(1,1),Sigma_error(2,1),Sigma_error(2,2)});
P=double(P);  

epshatsq=zeros(3,3,Teffective);
        for t=1:Teffective
            vecerror=reshape(epshat(t,:)'*epshat(t,:),[],1);
            vecerror(2)=[];

            temp=vecerror*vecerror';
       epshatsq(:,:,t)=temp;
        end
        
            vecsigma=reshape(Sigma_error,[],1);
            vecsigma(2)=[];
 Sigma_sigma=mean(epshatsq,3)-vecsigma*vecsigma';
  
  Beff=zeros(q,q,maxh);
varcovirf=zeros(q*q,q*q,maxh);
s_errors=zeros(q,q,maxh);

for h=1:maxh
      Hold=VAR_Companion^h;
Beff(:,:,h)=double(Hold(1:q,1:q)*P);

    G=Sigma_error;

    if h>1
for j=1:h-1
   G=G+(Beff(:,:,j)*Sigma_error*Beff(:,:,j)');
end
    end

    Hold=VAR_Companion^h;
Hold=Hold(1:q,1:q);

varcovirf(:,:,h)=(kron(eye(q),G)/Teffective)+(kron(eye(q),Hold)*H*Sigma_sigma*H'*kron(eye(q),Hold')/Teffective);

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