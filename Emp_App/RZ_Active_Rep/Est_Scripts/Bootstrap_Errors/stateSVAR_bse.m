function bs_irf_dist = stateSVAR_bse(y_tvar,I,beta_a,beta_b,p, hor,B, sind, rind)
% STATESVAR_BSE: computes Bootstrap Standard Errors for SVAR IRFs
% Inputs:
    % y_tvar is the vector of tvar observations (T x k)
    % I is the indicator variable (=1 for A, =0 other)
    % beta_a is the estimated beta for state A
    % beta_b is the estimated beta for state B
    % res is the matrix of size (T-p) x k
    % p is the VAR order
    % hor is the horizon for IRFs
    % B is number of bootstrap replications
    % sind is index of irf impulse variable in y
    % rind are indices of irf response variables in y
% Outputs:
    % bs_irf_dist is B x r x hor x 2 (where r is no. of IRF response variables)
    
    % setup
    [T,K] = size(y_tvar);
    TT = T-p; % truncated p is the length of fitted y and residuals
    r = length(rind);

    % Step 0: compute the residual matrix
    yhat_tvar = zeros(TT,K);
    X = lagmatrix(y_tvar, [1:p]);
    X = [ones(TT,1) X(p+1:end,:)];
    tI = I(p+1:end,:);

    for i=1:TT
        yhat_tvar(i,:) = (X(i,:)*tI(i))*beta_a + (X(i,:)*(1-tI(i)))*beta_b;
    end

    ty_tvar = y_tvar(p+1:end,:);
    resid = ty_tvar - yhat_tvar;
    resid = resid - mean(resid);

    % create bs placeholder matrix
    bs_irf_dist = zeros(B,r,hor,2);

    for b=1:B
        y_bs = zeros(T,K); % we compute a T x K dependent variable via RFVAR DGP 
        y_bs(1:p,:) = y_tvar(1:p,:);
        res_samp = datasample(resid, TT, 'Replace', true);
    
        for j=p+1:T
            ind = I(j);
            lagy_bs = [];
            for l=1:p
                lagy_bs = [lagy_bs, y_bs(j-l,:)];
            end
            x_bs = [1, lagy_bs];
            res = res_samp(j-p,:);
            y_bs(j,:) = (x_bs*ind)*beta_a + (x_bs*(1-ind))*beta_b + res;
        end

        % (estimate TVAR irfs using the SVAR func)
        [svara_temp, svarb_temp, ~, ~, ~, ~] = stateSVAR(y_bs,I,p,hor,rind,sind, B, 0, 0, 'normal');
        
        % (store in irf distribution)
        bs_irf_dist(b,:,:,1) = svara_temp;
        bs_irf_dist(b,:,:,2) = svarb_temp;
    end
end

