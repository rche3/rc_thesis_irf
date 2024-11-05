function [CoverageLPHARtemp, LengthLPHARtemp, BiasLPHARtemp] = lp_har_summary(y,arp,maxh,confint,trueirf) 

%Calculates coverage, length, and bias for LP OLS estimator


%Input:
    %y is the qxT matrix of variables
    %arp is the lag length
    %nstraps is the number of bootstrap replications
    %confit in the size of the desired confidence interval (e.g. .95)
    %trueirf are the actual Wold irfs from the data generating process

%Output:
    %CoverageLPHARtemp is a (q*q)x(maxh) matrix of coverage estimates (1 if confidence interval contains true estimate; 0 otherwise)
    %LengthLPHARtemp (q*q)x(maxh) matrix of length estimates
    %BiasLPHARtemp (q*q)x(maxh) matrix of bias estimates

Y=y';  [q T]=size(Y); 
saveY=Y; 

CoverageLPHARtemp=zeros(maxh,q*q);
LengthLPHARtemp=zeros(maxh,q*q);
BiasLPHARtemp=zeros(maxh,q*q);

alpha=1-confint;
upper=1-(alpha/2);

for h=1:maxh
  
    p=1+arp*q; 
F=zeros(q*arp,T-arp-h+1); 
for j=1:arp 
    F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-h+1-j)); 

end

saveF=F;

yeff=saveY(:,arp+h:end)'; %Dependent variable

Teffective=size(yeff',2);
F=[ones(1,Teffective);saveF]; %Design matrix

Hold_beta=zeros(q,q);
Hold_se=zeros(q,q);


for i=1:q
    
[betahat,~,se_beta,~,~,q_df] = har_ols(yeff(:,i),F',0);

Hold_beta(i,:)=betahat(2:q+1)';
Hold_se(i,:)=se_beta(2:q+1)';

end

    
 lpnw_intervals=[reshape(Hold_beta,[],1)-tinv(upper,q_df)*reshape(Hold_se,[],1) reshape(Hold_beta,[],1)+tinv(upper,q_df)*reshape(Hold_se,[],1)]; 
 
                  LengthLPHARtemp(h,:)=(lpnw_intervals(:,2)-lpnw_intervals(:,1))';

            if q==1
                                    jack= (reshape(trueirf(h),[],1) >= lpnw_intervals(:,1)) .* (reshape(trueirf(h),[],1) <= lpnw_intervals(:,2)); 
                                    
                                                                  CoverageLPHARtemp(h,(i-1)*q+(1:q))=jack';

                    
                             BiasLPHARtemp(h,(i-1)*q+(1:q))=betahat(2:q+1)-trueirf(h)'; 

            else
                    jack= (reshape(trueirf(:,:,h),[],1) >= lpnw_intervals(:,1)) .* (reshape(trueirf(:,:,h),[],1) <= lpnw_intervals(:,2)); 
                    
                              CoverageLPHARtemp(h,:)=jack';

                    
                             BiasLPHARtemp(h,:)=(reshape(Hold_beta,[],1)-reshape(trueirf(:,:,h),[],1))'; 

                    
            end
          
         
end
          

return