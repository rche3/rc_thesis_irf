function bs_irf_dist = linSVAR_bse(y,beta,resid,p,hor,B,sind,rind,normalise)
% LINSVAR_BSE: computes Bootstrap Standard Errors for Linear SVAR IRFs
% Inputs:
    % y is the vector of observations (T x k)
    % beta is the estimated beta
    % resid is the matrix of size (T-p) x k
    % p is the VAR order
    % hor is the horizon for IRFs
    % B is number of bootstrap replications
    % sind is index of irf impulse variable in y
    % rind are indices of irf response variables in y
% Outputs:
    % bs_irf_dist is B x r x hor (where r is no. of IRF response variables)
    
    % setup
    [T,K] = size(y);
    TT = T-p; % truncated p is the length of fitted y and residuals
    rsize = length(rind);

    % Step 0: centre the residuals
    resid = resid - mean(resid);

    % create bs placeholder matrix
    bs_irf_dist = zeros(B,rsize,hor);

    for b=1:B
        % disp(['Bootstrap rep: ', num2str(b), ' out of ', num2str(B)])
        y_bs = zeros(T,K); % we compute a T x K dependent variable via RFVAR DGP 
        y_bs(1:p,:) = y(1:p,:);
        res_samp = datasample(resid, TT, 'Replace', true);
        for j=p+1:T
            lagy_bs = [];
            for l=1:p
                lagy_bs = [lagy_bs, y_bs(j-l,:)];
            end
            x_bs = [1, lagy_bs];
            res = res_samp(j-p,:);
            y_bs(j,:) = x_bs*beta + res;
        end

        % (estimate TVAR irfs using the SVAR func)
        [bs_irf_temp, ~] = linSVAR(y_bs,hor,p,rind,sind,0,0,0,normalise); % clevel, bootstrap, B all set to zero as not used
        
        % (store in irf distribution)
        bs_irf_dist(b,:,:) = bs_irf_temp;
    end
end

