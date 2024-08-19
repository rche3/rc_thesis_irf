function [bs_beta_se, bs_beta_mean] = compute_block_bootstrap_se(resid, y, x, beta, h, nlag, y_pos_control, rpos)
% COMPUTE_BLOCK_BOOTSTRAP_SE Computes block bootstrapped standard errors for LP regression
% Returns the block bootstrapped standard errors for each beta in the LP
% regression given horizon, observations, and residuals

% Check if sizes are the same for resid, y, x
if size(resid, 1) == size(y, 1) && size(x, 1) == size(y,1)
    % pass
else
    error('Error: dimensions of resid, regressors, and dep. var not equal')
end

% Bootstrap settings
bs_samples = 1000;
block_length = ceil(size(resid, 1)^(1/3)); % Rule of thumb for block length

% Setup
T = length(resid);
beta_bs = nan(size(beta, 1), bs_samples);
Y = [y, x];

% Create block bootstrap samples
for j = 1:bs_samples
    % Generate random starting points for blocks
    start_points = randi(T - block_length + 1, ceil(T / block_length), 1);
    
    % Create indices for the block bootstrap sample
    indices = [];
    for i = 1:length(start_points)
        indices = [indices; (start_points(i):min(start_points(i)+block_length-1, T))'];
    end
    indices = indices(1:T); % Trim to exactly T observations
    
    % Create the block bootstrap sample
    resid_sample = Y(indices, :);
    y = resid_sample(:, 1);
    x = resid_sample(:, 2:end);
    
    % Estimate bootstrapped beta
    results = nwest(y, x, 0); % Note: x is already lagged
    beta_bs(:, j) = results.beta;
end

% Compute variance of the bootstrapped betas
cov_beta = cov(beta_bs');
bs_beta_var = diag(cov_beta);
bs_beta_se = sqrt(bs_beta_var);
bs_beta_mean = mean(beta_bs, 2); % Returns the mean of the betas across each row

end