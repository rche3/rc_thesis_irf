function Beff = lp_gls_analytical_simp(y, arp, maxh) 

%Analytical LP GLS Code from Lusompa (2023), modified for point estimates only without bias adjustment

%Input:
    %y is the qxT matrix of variables
    %arp is the lag length
    %maxh is the maximum desired LP horizon

%Output:
    %Beff is the qxqxmaxh matrix of the wold impulse responses. For
    %example Beff(:,:,h)' gives the horizon h estimate of the wold impulse

Y = y; 
[q, T] = size(Y); 
  
%population estimates for h=1
h = 1;
saveY = Y; 
  
p = 1 + arp*q; 
F = zeros(q*arp, T-arp-h+1); 
for j = 1:arp 
    F((j-1)*q+(1:q), :) = Y(:, (arp-j+1):(T-h+1-j));  
end
  
yeff = saveY(:, arp+h:end)'; %Dependent variable
 
Teffective = size(yeff', 2);
F = [ones(1, Teffective); F]; % Design Matrix
 
[betahat, ~, ~, epshat] = ols(yeff, F', 0);

Beff = zeros(q, q, maxh);
Beff(:, :, h) = betahat(2:q+1, :);

epshat = epshat - mean(epshat);

%population estimate for h>=2 
for h = 2:maxh 
    epscorrect = zeros(T-arp-h+1, q);  
          
    %gls correction for population estimate
    for hsub = 1:h-1 
        epscorrect = epscorrect + epshat(h-hsub:T-arp-hsub, :) * Beff(:, :, hsub);
    end 
    
    Y = y;  
    saveY = Y; 
 
    yeff = saveY(:, arp+h:end)' - epscorrect; %Dependent Variable

    p = 1 + arp*q; 
    F = zeros(q*arp, T-arp-h+1); 
    for j = 1:arp 
        F((j-1)*q+(1:q), :) = Y(:, (arp-j+1):(T-h+1-j)); 
    end
  
    Teffective = size(yeff', 2);
    F = [ones(1, Teffective); F]; % Design Matrix
 
    %population estimate
    [betahat, ~, ~, ~] = ols(yeff, F', 0);
    
    Beff(:, :, h) = betahat(2:q+1, :);
end

return