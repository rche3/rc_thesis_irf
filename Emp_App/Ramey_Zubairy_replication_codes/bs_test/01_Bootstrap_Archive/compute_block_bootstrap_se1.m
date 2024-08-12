function [bs_beta_se, bs_beta_mean] = compute_block_bootstrap_se1(resid, y, x, beta)
% COMPUTE_BLOCK_BOOTSTRAP_SE Computes block bootstrapped standard errors for LP regression
% Returns the block bootstrapped standard errors for each beta in the LP
% regression by resampling blocks of residuals

% Check if sizes are the same for resid, y, x
if size(resid, 1) == size(y, 1) && size(x, 1) == size(y,1)
    % pass
else
    error('Error: dimensions of resid, regressors, and dep. var not equal')
end

% Bootstrap settings
bs_samples = 1000;
T = length(resid);
block_length = ceil(T^(1/3)); % Rule of thumb for block length

% Setup
demeaned_resid = resid - mean(resid);
beta_bs = nan(size(beta, 1), bs_samples);

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
    
    % Create the block bootstrap sample of residuals
    resid_sample = demeaned_resid(indices);
    
    % Create the bootstrap dependent variable
    y_bs = x * beta + resid_sample;
    
    % Estimate bootstrapped beta
    results = nwest(y_bs, x, 0); % Note: x is already lagged
    beta_bs(:, j) = results.beta;
end

% Compute standard errors of the bootstrapped betas
bs_beta_se = std(beta_bs, 0, 2);
bs_beta_mean = mean(beta_bs, 2); % Returns the mean of the betas across each row