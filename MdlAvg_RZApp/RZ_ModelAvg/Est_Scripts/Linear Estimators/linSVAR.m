function [liny,confidencey, beta, mult_lpu] = linSVAR(y,hor,p,rind,sind,clevel,bootstrap,B,normalise, method)
    % SVAR computes IRFs (Cholesky) and bootstrapped standard errors using SVAR estimation
    % Inputs:
        % y is the T x k vector of observations
        % hor is the number of IRF horizons
        % rind are indices of irf response variables
        % sind is the index of irf shock variables
        % clevel is the critical value for confidence intervals (default 95%)
        % bootstrap is whether to use bootstrap (i.e., do you want CIs?)
        % p is the number of lags
    % Outputs
        % liny is the r x hor matrix of irf responses
        % confidencey is the 2 x r x h matrix of upper and lower CI bounds
        % mult_lpu is as 3 x hor matrix of multipliers of liny(2,:) on liny(1,:)

    % setup
    [T,k] = size(y);
    rsize = length(rind);

    % create lagged values and drop obs
    ylag = lagmatrix(y, [1:p]);
    ylag = ylag(p+1:end,:);
    y = y(p+1:end,:);

    % create placeholder matrices
    liny = nan(rsize,hor);
    confidencey = nan(2,rsize,hor);
    mult_lpu = zeros(3,hor);

    % COMPUTE IRFS
    % Step 1: Completed reduced form regression
    TT = size(y,1);
    X = [ones(TT,1) ylag];
    beta = inv(X'*X)*X'*y;

    % Step 2: Obtain residuals and perform Cholesky decomposition
    resid = y - X * beta;
    sigma_u = (resid'* resid)/(TT-size(X,2));
    P = chol(sigma_u, 'Lower'); % P is equal to B^{-1}
%     disp(P)

    % Step 3: Compute the IRFs recursively
    theta = zeros(k,k,hor); % final k x k x hor matrix of impulse responses
    phi = zeros(k,k,hor); % VMA coefficients placeholder
    phi(:,:,1) = eye(k); % initialise phi contemporaneous with I
    
    % reshape the RF-VAR beta into accessible coeffs
    A = zeros(k,k,p);
    A_matrix = beta(2:end,:);
    for j=1:p
        Aj_ind = 1+(j-1)*k:j*k;
        A(:,:,j) = A_matrix(Aj_ind,:)'; % transpose since the beta matrix holds transposed A_1, .. A_p
    end
    
    % recursively compute orthogonalised IRFs
    for h=2:hor
        phi_temp = zeros(k,k);
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
    theta_i = zeros(rsize,hor); 
    for j=1:rsize
        resp_ind = rind(j);
        for h=1:hor
            theta_i(j,h) = theta(resp_ind,sind,h); % third arg is the shock index
        end
    end
    
    norm = theta(1,1,1); % normalise with respect to shock

    if normalise == 1
        theta_i = theta_i / norm;
    else
        % pass
    end
    liny = theta_i;
    
    if size(liny,1) > 1
        mult_lpu(3,:) = cumsum(liny(2,:))./cumsum(liny(1,:));
    else
    end

    % Step 5: Bootstrap inference
    if bootstrap == 1
        bs_irf_dist = linSVAR_bse(y,beta,resid,p,hor,B,sind,rind,normalise);
        switch method
            case "normal"
                se = zeros(rsize, hor);
                for ri=1:rsize
                    for hi=1:hor
                        bsvec = squeeze(bs_irf_dist(:,ri,hi));
        %                 disp([mean(bsvec), liny(ri,hi)])
                        se(ri,hi) = std(bsvec);
                    end
                end
                confidencey(1,:,:) = liny - clevel * se;
                confidencey(2,:,:) = liny + clevel * se;
            case "percentile"
                alpha = 1 - normcdf(clevel);
                for ri=1:rsize
                    for hi=1:hor
                        bsvec = squeeze(bs_irf_dist(:,ri,hi));
                        upperq = quantile(bsvec, 1-alpha);
                        lowerq = quantile(bsvec, alpha);
                        confidencey(1,ri,hi) = lowerq;
                        confidencey(2,ri,hi) = upperq;
                    end
                end
            otherwise 
                error('Invalid method specified. Use "normal" or "percentile"')
        end

        if size(liny,1) > 1
            mult_bs_dist = cumsum(squeeze(bs_irf_dist(:,2,:)),2)./cumsum(squeeze(bs_irf_dist(:,1,:)),2);
            switch method
                case "normal"
                    % multiplier CI
                    mult_se = var(mult_bs_dist,1,1).^(1/2);
                    mult_lpu(1,:) = mult_lpu(3,:) - clevel * mult_se;
                    mult_lpu(2,:) = mult_lpu(3,:) + clevel * mult_se;
                case "percentile"
                    % multiplier CI
                    alpha = 1-normcdf(clevel);
                    mult_lpu(1,:) = quantile(mult_bs_dist, 1-alpha);
                    mult_lpu(2,:) = quantile(mult_bs_dist, alpha);
            otherwise 
                error('Invalid method specified. Use "normal" or "percentile"')
            end
        else
        end
    else
        confidencey = zeros(2,rsize,hor);
    end

    % Step 6: Bootstrap inference for multiplier
    
end

