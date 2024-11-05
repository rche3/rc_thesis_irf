function [B, A, mu, opt_lag] = generate_DGP_rz(max_lag, difference)
    % Setup remains the same
    timesample = 1;
    statechoice = 1;
    shockchoice = 1;
    transformation = 1;
    taxchoice = 0;
    datafirstdiff = 0;
    trend = 0;
    opt = 0;
    
    % Initialize variables for lag selection
    aic = zeros(max_lag, 1);
    bic = zeros(max_lag, 1);

    % Compute AIC and BIC for different lag lengths
    for p = 1:max_lag
        nlag = p;

        run('figure_irfs_multipliers_setup.m')

        if difference == 1
                y = diff(y_tvar);
        else
            y = y_tvar;
        end

        [T, k] = size(y);
        ylag = lagmatrix(y, [1:p]);
        y_temp = y(p+1:end,:);
        ylag = ylag(p+1:end,:);
        [Teff, ~] = size(y_temp);
        
        % Create X matrix
        X = [ones(Teff, 1) ylag];
        
        % Estimate
        beta = X \ y_temp;
        
        % Compute residuals and covariance matrix
        resid = y_temp - X*beta;
        sigma = (resid'* resid)/Teff;
        
        % Compute information criteria
        n_params = k * (k*p + 1);  % Number of parameters
        log_det_sigma = log(det(sigma));
        aic(p) = log_det_sigma + (2/Teff)*n_params;
        bic(p) = log_det_sigma + (log(Teff)/Teff)*n_params;
    end
    
    % Select optimal lag based on AIC
    [~, opt_lag_aic] = min(aic);
    % Select optimal lag based on BIC
    [~, opt_lag_bic] = min(bic);
    
    % Use BIC optimal lag (tends to be more conservative)
    opt_lag = opt_lag_bic;
    
    % Estimate VAR with optimal lag
    ylag = lagmatrix(y, [1:opt_lag]);
    y = y(opt_lag+1:end,:);
    ylag = ylag(opt_lag+1:end,:);
    [Teff, k] = size(y);
    
    % Create X
    X = [ones(Teff, 1) ylag];
    
    % Estimate
    beta = X \ y;
    
    % Grab the residuals and compute B^-1
    resid = y - X*beta;
    sigma_u = (resid'* resid)/(Teff-size(X,2));
    Binv = chol(sigma_u, 'Lower'); % P is equal to B^{-1}
    B = inv(Binv);
    
    % Return the final values of the matrices
    mu = beta(1,:)';
    A = cell(1,opt_lag); % Preallocate cell array
    for i = 1:opt_lag
        A{i} = beta(2+(i-1)*k:1+i*k,:)'; % transpose to turn into A1, A2 etc. in normal math form
    end

    % (rescale B to have unit diagonals)
    B = B./diag(B);
    
    % Display lag selection results
    fprintf('Optimal lag length (AIC): %d\n', opt_lag_aic);
    fprintf('Optimal lag length (BIC): %d\n', opt_lag_bic);
    fprintf('Selected lag length: %d\n', opt_lag);
end