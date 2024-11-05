function Beff = lp_gls_analytical_simp(Y, nlag, maxh, recurShock, responseV) 

%Analytical LP GLS Code from Lusompa (2023), modified for point estimates only without bias adjustment

%Input:
    %y is the qxT matrix of variables
    %nlag is the lag length
    %maxh is the maximum desired LP horizon
    %recurShock is the location of the shock in the data vector
    %responseV is the location of the responseV in the data vector

%Output:
    %Beff is the (maxh)x1 vector of the specific shock-response pair IRFs
     % (note, usually maxh is already the actual irf horizon + 1 to account for matlab indexing from 1)

[Q, T] = size(Y);
Beff = zeros(maxh,1);
Yt = Y';

% Estimation for h=0 (instantaneous response, matches OLS LP)
h = 0;
% disp('This is the y data for LP_gen_data in my GLS implementation')
% disp(size(Yt))
% disp(Yt(1:10,:))
[y, x, w] = LP_gen_data(Yt, recurShock, responseV, nlag, h);  % Note: y is transposed here
% disp(size(y));
% disp(y(1:10,:));
X = [ones(size(x,1),1), x, w];
[Beta, ~, ~, ~] = LS(y, X);
Beff(h+1) = Beta(2);  % The impulse response coefficient

% Estimation for h>=1
for h = 2:maxh-1
    F = zeros(T-nlag-h+1, Q*nlag);
    y_temp = Yt;
    for j = 1:nlag
        F(:,(j-1)*Q+(1:Q)) = y_temp((nlag-j+1):(T-h+1-j), :);
    end
    
    yeff = y_temp((nlag+h):end, responseV);
    Teffective = size(yeff, 1);
    F = [ones(Teffective, 1), F];
    
    epscorrect = zeros(T-nlag-h+1, 1);
    for hsub = 1:h-1
        epscorrect = epscorrect + y_temp((nlag+h-hsub):(T-hsub), recurShock) * Beff(hsub+1);
    end
    yeff = yeff - epscorrect;
    
    [betahat, ~, ~, ~] = LS(yeff, F);
    Beff(h+1) = betahat(recurShock+1);
end

end