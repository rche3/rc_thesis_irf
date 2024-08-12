function [bs_beta_se, bs_beta_mean] = compute_bootstrap_se(resid, y, x, beta, h, nlag, y_pos_control, wild)
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
bs_samples = 200;

% setup
T_resid = length(resid); % by design, this should be T - h + 1 in length, h \in [1,20]
demeaned_resid = resid - sum(resid)/T_resid;
beta_bs = nan(size(beta, 1), bs_samples);
T = T_resid - h - nlag; % our effective "T" which will be the length of the bootstrapped dependent variable
y_bs = zeros(T, 1);

% wild bootstrap settings
e1 = (-sqrt(5)-1)/2;
e2 = (sqrt(5)+1)/2;
p1 = (sqrt(5)+1)/(2*sqrt(5));
p2 = 1-p1;
values = [e1, e2];
probabilities = [p1, p2];
num_samples = T;

% create residual bootstrap samples

for j = 1:bs_samples
% Print to the Command Window
    resid_sample = datasample(demeaned_resid, T, 'Replace', true);

    % create the bootstrap dependent variable
    x_temp = x; % this is a temporary x variable which will become iteratively updated, initialised as our regressor matrix

    % create the wild bootstrap errors
    wild_errors = disc_rv(values, probabilities, num_samples); % randomly sampled from discrete pd

    % DWB residuals / disturbances
    

    for i = 1:T
        x_temp = x;
        
        % compute and store y_bs variables
        if wild == 1
            resid = resid_sample(i) * wild_errors(i);
        else
            resid = resid_sample(i);
        end
        y_bs(i) = x_temp(i, :) * beta + resid; 
        controls = x_temp(:, 3:end);
        controls(i, y_pos_control) = y_bs(i); % each y_bs should override the control in x_temp
        
        if i > 1
            for n = 2:min(i, 4)
                % ensure each 1x3 control is one version lagged of the prior
                lag = n-1;
                read_pos = 3 * (lag-1) + y_pos_control;
                write_pos = 3 * lag + y_pos_control;
                controls(i-lag, write_pos) = controls(i, read_pos);
            end
        end
        x_temp = [x_temp(:, 1:2), controls];
%             from i = nlags +  h + 1 onwards, the entire x_temp control vector
%             as the regressor should be bootstrap variables
    end

    % estimate bootstrapped beta
    x_trunc = x_temp(1:T, :);
    results=nwest_rc(y_bs, x_trunc, 0); % note we don't need to lag x as input arg x is already lagged
    beta_bs(:, j) = results.beta;
end

% compute variance of the bootstrapped betas
cov_beta = cov(beta_bs');
bs_beta_var = diag(cov_beta);
bs_beta_se = bs_beta_var.^(1/2);

bs_beta_mean = mean(beta_bs, 2); % returns the mean of the betas across each row
% fprintf('bs beta mean for horizon %d is %.4f \n', h, bs_beta_mean(rpos));

end

