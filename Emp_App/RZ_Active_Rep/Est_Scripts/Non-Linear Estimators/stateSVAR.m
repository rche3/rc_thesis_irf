function [irfa, irfb, cia, cib, bs_beta_dist, pval, beta] = stateSVAR(y, I, p, hor, rind, sind, B, bootstrap, clevel, method)
% STATESVAR computes the irfs and CIs for the state dependent vector autoregression
% By default uses the Cholesky decomposition for SVAR identification
% Inputs:
    % y is the complete vector of observations T x k
    % I is a T x 1 indicator variable (taking the value 1 for State "A" and 0 for State "B")
    % p is the order of the VAR
    % hor is the IRF horizon
    % rind are the indices of the variables in y to compute IRFs for
    % sind is the shock location in k-dimensional time series
    % B is the number of bootstrap reps
    % bootstrap is 1 or 0 (compute bootstrap SEs or just leave blank to prevent infinite recursion)
    % clevel is the confidence level e.g. 1.96 for 5% two sided
    % method is the method for calculating CIs
% Outputs:
    % irfa is the (2 x hor) impulse response for State A
    % irfb is the (2 x hor) impulse response for State B
    % cia is the CIfor State A
    % cib is the CI for State B
    % pvala is p-value for the State A estimated IRFs (against null of zero)
    % pvala is p-value for the State B estimated IRFs (against null of zero)
    % beta is a 2 x (Kp+1) x K matrix of residuals

    %%%% FUNCTION START %%%%
    
    % SETUP
    % create lagged variables and drop observations
    ylag = lagmatrix(y, [1:p]);
    ylag = ylag(p+1:end,:);
    y = y(p+1:end,:);
    I = I(p+1:end,:);
    
    % create the State A and B dependent vars
    ya = y(I==1,:); xa = ylag(I==1,:);
    yb = y(I==0,:); xb = ylag(I==0,:);
    [~, nvara] = size(ya);
    K = nvara;
    
    % create placeholder matrices
    states = {'a','b'}; num_states = length(states);  % assuming 'state' is already defined
    beta = zeros(K*p+1, K, num_states);
    resid = cell(1,num_states);
    sigma_u = cell(1,num_states);

    % output placeholders
    r = length(rind);
    irf_comb = zeros(length(states),r,hor);
    pval = nan(2,r,hor);
    bs_beta_dist = nan(2,B,r,hor);
    cia = nan(2,r,hor);
    cib = nan(2,r,hor);

    % COMPUTE IRFS
    for i=1:2
        % Step 1: Complete the reduced form regression
        yi = ['y' states{i}]; Y = eval(yi);
        xi = ['x' states{i}]; XX = eval(xi);
        TT = size(Y, 1);
        X = [ones(TT,1) XX];
        beta_temp = inv(X'*X) * X' * Y; % (Kp+1) x K

        % Step 2: Obtain residuals and perform Cholesky decomposition
        resid_temp = Y - X * beta_temp;
        sigma_u_temp = (resid_temp'* resid_temp)/(TT-size(XX,2));
        P = chol(sigma_u_temp, 'Lower'); % P is equal to B^{-1}

        % Step 3: Compute the IRFs recursively
        theta = zeros(K,K,hor); % final KxKxH matrix of impulse responses
        phi = zeros(K,K,hor); % VMA coefficients placeholder
        phi(:,:,1) = eye(K); % initialise phi contemporaneous with I
        
        % reshape the RF-VAR beta into accessible coeffs
        A = zeros(K,K,p);
        A_matrix = beta_temp(2:end,:);
        for j=1:p
            Aj_ind = 1+(j-1)*K:j*K;
            A(:,:,j) = A_matrix(Aj_ind,:)'; % transpose since the beta matrix holds transposed A_1, .. A_p
        end
        % recursively compute orthogonalised IRFs
        for h=2:hor
            phi_temp = zeros(K,K);
            for s = 1:min(h-1,p)
                phi_temp = phi_temp + A(:,:,s)*phi(:,:,h-s);
            end            
            phi(:,:,h) = phi_temp;
        end

        for h=1:hor
            phi_slice_h = squeeze(phi(:,:,h));
            theta(:,:,h) = phi(:,:,h) * P;
        end

        % Step 4: Pull specific IRFs to response vars & normalise
        theta_i = zeros(r,hor); 
        for j=1:r
            resp_ind = rind(j);
            for h=1:hor
                theta_i(j,h) = theta(resp_ind,sind,h); % third arg is the shock index
            end
        end

        norm = theta(1,1,1); % normalise with respect to shock
        theta_i = theta_i / norm;
        irf_comb(i,:,:) = theta_i;

        % (optional) store auxiliary variables
        beta(:,:,i) = beta_temp;
    end

    % retrieve IRFs
    irfa = squeeze(irf_comb(1,:,:));
    irfb = squeeze(irf_comb(2,:,:));

    % Step 5: Bootstrap inference
    if bootstrap == 1
        beta_a = beta(:,:,1); beta_b = beta(:,:,2);
        irf_bs_beta = stateSVAR_bse(y, I, beta_a, beta_b, p, hor, B, sind, rind);
        
        bsirf_a = squeeze(irf_bs_beta(:,:,:,1));
        bsirf_b = squeeze(irf_bs_beta(:,:,:,2));
        
        switch method
            case "normal"
                sea = zeros(r,hor);
                seb = zeros(r,hor);
                for ri=1:r
                    for hi=1:hor
                        bsvec_a = squeeze(bsirf_a(:,ri,hi));
                        sea(ri,hi) = std(bsvec_a);
                        bsvec_b = squeeze(bsirf_b(:,ri,hi));
                        seb(ri,hi) = std(bsvec_b);
                    end
                end
                cia(1,:,:) = irfa - clevel * sea;
                cia(2,:,:) = irfa + clevel * sea;
                cib(1,:,:) = irfb - clevel * seb;
                cib(2,:,:) = irfb + clevel * seb;
            case "percentile"
                % Calculate cia and cib for percentile method
                alpha = 1 - normcdf(clevel);
                for ri=1:r
                    for hi=1:hor
                        bsvec_a = squeeze(bsirf_a(:,ri,hi));
                        bsvec_b = squeeze(bsirf_b(:,ri,hi));
                        cia(1,ri,hi) = quantile(bsvec_a, alpha);
                        cia(2,ri,hi) = quantile(bsvec_a, 1-alpha);
                        cib(1,ri,hi) = quantile(bsvec_b, alpha);
                        cib(2,ri,hi) = quantile(bsvec_b, 1-alpha);
                    end
                end
            otherwise
                error('Invalid method specified. Use "normal" or "percentile"')
        end        
    else
        cia(1,:,:) = nan(size(cia,2), size(cia,3));
        cia(2,:,:) = nan(size(cia,2), size(cia,3));
        cib(1,:,:) = nan(size(cib,2), size(cib,3));
        cib(2,:,:) = nan(size(cib,2), size(cib,3));
    end
end

