function bs_beta_se = linlplagaug_dwbse(resid, y, x, beta)
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
resid = resid - sum(resid)/T;
beta_bs = nan(size(beta, 1), bs_samples);

% create residual bootstrap samples
for j = 1:bs_samples
    resid_sample = datasample(resid, T, 'Replace', true);

    % create the bootstrap dependent variable
    y_bs = x * beta + resid_sample;

    % estimate bootstrapped beta
    results=nwest(y_bs, x, 0); % note we don't need to lag x as input arg x is already lagged
    beta_bs(:, j) = results.beta;
end

% compute variance of the bootstrapped betas
% cov_beta = cov(beta_bs');
% bs_beta_vars = diag(cov_beta);
% bs_beta_se = bs_beta_vars.^(1/2);

bs_beta_se = std(beta_bs, 0, 2);

end

