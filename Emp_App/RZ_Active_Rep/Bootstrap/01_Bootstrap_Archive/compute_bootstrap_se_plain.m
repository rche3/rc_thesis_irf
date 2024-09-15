function [bs_beta_se, bs_beta_mean] = compute_bootstrap_se_plain(resid, y, x, beta, h, nlag, y_pos_control, rpos)
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
bs_samples = 1000;

% setup
T = length(resid); % by design, this should be T - h + 1 in length, h \in [1,20]
beta_bs = nan(size(beta, 1), bs_samples);
Y = [y, x];

% create residual bootstrap samples
for j = 1:bs_samples

    resid_sample = datasample(Y, T, 'Replace', true);
    y = resid_sample(:, 1);
    x = resid_sample(:, 2:end);

    % estimate bootstrapped beta
    results=nwest(y, x, 0); % note we don't need to lag x as input arg x is already lagged
    beta_bs(:, j) = results.beta;
end

% compute variance of the bootstrapped betas
cov_beta = cov(beta_bs');
bs_beta_var = diag(cov_beta);
bs_beta_se = bs_beta_var.^(1/2);

bs_beta_mean = mean(beta_bs, 2); % returns the mean of the betas across each row
% fprintf('bs beta mean for horizon %d is %.4f \n', h, bs_beta_mean(rpos));

end

