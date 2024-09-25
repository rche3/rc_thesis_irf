function [liny,confidencey] = linSVAR(y,hor,p,rind,sind,clevel,bootstrap,B,normalise)
    % SVAR computes IRFs (Cholesky) and bootstrapped standard errors using SVAR estimation
    % Inputs:
        % y is the T x k vector of observations
        % hor is the number of IRF horizons
        % rind are indices of irf response variables
        % sind is the index of irf shock variables
        % clevel is the critical value for confidence intervals (default 95%)
        % p is the number of lags
    % Outputs
        % liny is the r x hor matrix of irf responses
        % confidencey is the 2 x r x h matrix of upper and lower CI bounds

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

    % COMPUTE IRFS
    % Step 1: Completed reduced form regression
    TT = size(y,1);
    X = [ones(TT,1) ylag];
    beta = inv(X'*X)*X'*y;

    % STep 2: Obtain residuals and perform Cholesky decomposition
    resid = y - X * beta;
    sigma_u = (resid'* resid)/(TT-size(X,2));
    P = chol(sigma_u, 'Lower'); % P is equal to B^{-1}

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

    % Step 5: Bootstrap inference
    if bootstrap == 1
        irf_bs_beta = linSVAR_bse(y,beta,resid,p,hor,B,sind,rind,normalise);
        se = zeros(rsize, hor);
        for ri=1:rsize
            for hi=1:hor
                bsvec = squeeze(irf_bs_beta(:,ri,hi));
                se(ri,hi) = std(bsvec);
            end
        end
        confidencey(1,:,:) = liny - clevel * se;
        confidencey(2,:,:) = liny + clevel * se;
    else
        confidencey = zeros(2,rsize,hor);
    end
end

