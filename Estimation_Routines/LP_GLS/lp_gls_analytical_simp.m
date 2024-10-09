function Beff = lp_gls_analytical_simp(y, arp, maxh, recurShock, responseV) 

%Analytical LP GLS Code from Lusompa (2023), modified for point estimates only without bias adjustment

%Input:
    %y is the qxT matrix of variables
    %arp is the lag length
    %maxh is the maximum desired LP horizon
    %recurShock is the location of the shock in the data vector
    %responseV is the location of the responseV in the data vector

%Output:
    %Beff is the (maxh+1)x1 vector of the specific shock-response pair IRFs

[Q, T] = size(y);
Beff = zeros(maxh+1,1);

% Estimation for h=1 (instantaneous response, matches OLS LP)
h = 1;
[Y, X, W] = LP_gen_data(y', recurShock, responseV, arp, h);  % Note: y is transposed here
X_full = [ones(size(X,1),1), X, W];
[Beta, ~, ~, ~] = LS(Y, X_full);
Beff(h) = Beta(2);  % The impulse response coefficient

% Estimation for h>=1
for h = 2:maxh
    F = zeros(T-arp-h+1, Q*arp);
    y_temp = y';
    for j = 1:arp
        F(:,(j-1)*Q+(1:Q)) = y_temp((arp-j+1):(T-h+1-j), :);
    end
    
    yeff = y_temp((arp+h):end, responseV);
    Teffective = size(yeff, 1);
    F = [ones(Teffective, 1), F];
    
    epscorrect = zeros(T-arp-h+1, 1);
    for hsub = 1:h-1
        epscorrect = epscorrect + y_temp((arp+h-hsub):(T-hsub), recurShock) * Beff(hsub+1);
    end
    yeff = yeff - epscorrect;
    
    [betahat, ~, ~, ~] = LS(yeff, F);
    Beff(h+1) = betahat(recurShock+1);
end

end