function Beff = lp_gls_analytical_manual(Y, nlag, maxh, recurShock, responseV, fgls_method) 

% Analytical LP GLS Code using the method described in Lusompa (2023)

%Input:
    %y is the Txq matrix of variables
    %nlag is the lag length
    %maxh is the maximum desired LP horizon
    %recurShock is the location of the shock in the data vector
    %responseV is the location of the responseV in the data vector

%Output:
    %Beff is the (maxh)x1 vector of the specific shock-response pair IRFs
     % (note, usually maxh is already the actual irf horizon + 1 to account for matlab indexing from 1)

[T, Q] = size(Y);
Beff = zeros(maxh,1);

% Estimation for h=0 (instantaneous response, matches OLS LP)
h=0;
[y, x, w] = LP_gen_data(Y, recurShock, responseV, nlag, h);
X = [ones(size(x,1),1), x, w];
[Beta, ~, ~, ~] = LS(y, X);
Beff(h+1) = Beta(2);  % The impulse response coefficient

% Estimation for h >= 2 (from the 3rd element in the Beff vector)
% (first we need to estimate the regular OLS residuals for h=1 to h=H-1, epsilon_{t+i,k} is Lutkepohl paper)
% eps = nan(T,maxh-1); % create a (maxh-1) x 1 empty matrix for storage of epsilon_{t+1,k} to epsilon_{t+H,k}
% beta = nan(maxh-1,1);
% for h = 1:maxh-1
%     [y, x, w] = LP_gen_data(Y, recurShock, responseV, nlag, h);
%     X = [ones(size(x,1),1), x, w];
%     [betahat, ~, ~, res] = LS(y, X);
%     eps(:,h) = [res; nan(nlag+h, 1)];
%     beta(h) = betahat;
% end

% to perform all later FGLS transformations to y, we just use epshat, which is \epsilon_{t+1,k}
if maxh > 2 % if the max mathematical horizon is H = 1
    % first compute the first horizon OLS regression to obtain epsilon
    h = 1;
    [y1, x1, w1] = LP_gen_data(Y, recurShock, responseV, nlag, h);
    X1 = [ones(size(x1,1),1), x1, w1];
    [Beta, ~, ~, Res] = LS(y1, X1);
    Beff(h+1) = Beta(2);  % The impulse response coefficient
    epshat = Res;
%     disp(size(epshat))
    [Teps, ~] = size(epshat);

    for h = 2:maxh-1
        % generate data first
        H = maxh-1;
        [y, x, w] = LP_gen_data(Y, recurShock, responseV, nlag, h);
        [Teff,~] = size(y);

        if fgls_method == 1 % Lusompa method
            epscorrect = zeros(Teff, 1);
            for hsub = 1:h-1
                eps_temp = epshat(Teps-(h-hsub)-Teff+1:Teps-(h-hsub),:); % this will be the Teffx1 vector of epsilon_{t+hsub}
                epscorrect = epscorrect + Beff(h-hsub+1) * eps_temp;
            end
        
            % generate the "y" and regressors as normal and then correct to \tilde{Y}
            X = [ones(size(x,1),1), x, w];
            ytilde = y - epscorrect;
            [betahat, ~, ~, ~] = LS(ytilde, X);
            Beff(h+1) = betahat(recurShock+1); % save the position of the coefficient estimate
        else % otherwise, it will be the Breitung method
            ytilde = y - epshat(Teps - Teff +1:Teps); % grab the final Teff observations of the eps_hat
            eps_regressors = zeros(Teff, h-2);
            for i=2:h-1
                delta_th = h-i;
                eps_regressors(:,i-1) = epshat(Teps-delta_th-Teff+1:Teps-delta_th);
            end
            X = [ones(size(x,1), 1), x, w, eps_regressors];
%             disp(eps_regressors)
            [betahat, ~, ~, ~] = LS(ytilde, X);
            Beff(h+1) = betahat(recurShock+1); % save the position of the coefficient estimate
        end
    end
end
end