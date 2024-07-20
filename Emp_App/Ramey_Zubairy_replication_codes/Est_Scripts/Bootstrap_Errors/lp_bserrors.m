function bs_beta_se = lp_bserrors(resid, y, x, beta)
% COMPUTE_BOOTSTRAP_SE Summary of this function goes here
% Returns the bootstrapped standard errors for each beta in the LP
% regression given horizon, observations, Newey-West residuals

% check the size is the same for resid, y, x
if size(resid, 1) == size(y, 1) && size(x, 1) == size(y,1)
    % pass
else
    disp('Error, dimensions of resid, regressors, and dep. var not equal')
end

% bootstrap settings
run('bootstrap_settings.m')

% setup
T = length(resid);
resid_sample = nan(T, bs_samples);
beta_bs = nan(size(beta, 1), bs_samples);

% create residual bootstrap samples
for i = 1:bs_samples
    resid_sample(:, i) = datasample(resid, T, 'Replace', true);

    % create the bootstrap dependent variable
    y_bs = nan(T, 1);
    for j = 1:T
        y_bs(j) = x(j, :) * beta + resid_sample(j, i); 
    end

    % estimate bootstrapped beta
    results=nwest(y_bs, x,i); % note we don't need to lag x as input arg x is already lagged
    beta_bs(:, i) = results.beta;
end

% compute variance of the bootstrapped betas
cov_beta = cov(beta_bs');
bs_beta_vars = diag(cov_beta);
bs_beta_se = bs_beta_vars.^(1/2);

end

