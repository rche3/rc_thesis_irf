function [Structural_IRF,fstat] = lp_gls_boot_iv(y,arp,nstraps,maxh,blocksize,badj) 

%LP GLS IV Bootstrap Code

%Input:
    %y is the qxT matrix of variables; the first row should be the
        %instrument and the second variable should be the variable that's
        %being instrumented
    %arp is the lag length
    %nstraps is the number of bootstrap replications
    %maxh is the maximum desired LP horizon
    %blocksize is the blocksize of the bootstrap
    %badj equals 1 for bias adjustment and 0 otherwise

%Output:
    %Structural_IRF qx(maxh) by nstraps matrix structural impulse responses
    %fstat is the bootstrap first stage F-statistic.

[Beff,Sigma_hold] = lp_gls_boot(y,arp,nstraps,maxh,blocksize,badj); 

q=size(y,1);


  %%%%%%%%%%%%%Structural Impulse Responses Calculation

  %LDL decomposition
Recursive_decomp=zeros(q,q,nstraps);
for g=1:nstraps
    L_chol=chol(Sigma_hold(:,:,g),'lower');
    diag_sqrt_L=(diag(diag(L_chol)));
    Recursive_decomp(:,:,g)=L_chol/(diag_sqrt_L); %L in LDL Decomposition

end


Structural_IRF=zeros(q,maxh+1,nstraps);
Structural_IRF_0=squeeze(Recursive_decomp(:,1,:)); %impact irf column of instrument
Structural_IRF(:,1,:)=Structural_IRF_0;

%impulse responses (for all variables at all horizons) to a change in the instrument
for horizon=1:maxh

   for g=1:nstraps
              
    Hold=Beff(:,:,horizon,g)';
           
    Structural_IRF_temp=Hold*Recursive_decomp(:,:,g); 
    Structural_IRF(:,horizon+1,g)=Structural_IRF_temp(:,1); 
   end  
end


%impulse responses (for all variables at all horizons) to a change in the
%instrument divided by the response of the independent variable of interest
%responses to the instrument; e.g. in the empirical application this is the
%responses of all variables at all horizons divided by the cotemporaneous
%responses of the one-year treasury
for horizon=1:maxh+1

   for g=1:nstraps
    Structural_IRF(:,horizon,g)=Structural_IRF(:,horizon,g)./Structural_IRF_0(2,g); 
   end  
end
 
%f stat calculation
                tstat=mean(Structural_IRF_0(2,:))/std(Structural_IRF_0(2,:));
                fstat=tstat.^2;

return
