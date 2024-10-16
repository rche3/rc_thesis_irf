function [IRF, n_lags_est] = resid_est(data_sim,settings)
%RESID_EST Summary of this function goes here
%   Detailed explanation goes here

run('Estimation_Setup'); % common setup for all estimation methods

%%% compute the LP residuals (LP errors) at max horizon H
res = IRF_LP_res(Y,recursiveShock,nlags,IRF_hor - 1); % residual matrix for LP-form of VAR(p)

theta_1_tH = res(:, policyV-1); % this is theta_1 per BL 2022
theta_i_tH = res(:, responseV-1); % this is theta_i, residuals for the response variable

%%% compute the estimator using the proxy z_t
% restrict the z_t to p+1:T-H
z_t = Y(nlags+1:end,1)';

window_size = IRF_hor;
num_windows = length(z_t) - window_size +1;

z_t_matrix = fliplr(hankel(z_t(1:num_windows), z_t(num_windows:end)));

% multiply each row of z_t_matrix by the corresponding scalar resid
numerator = sum(bsxfun(@times, theta_i_tH, z_t_matrix), 1);
% ^bsxfun is binary singleton expansion function - applies
% element-by-element operation to two arrays, in this case, for each "t",
% multiplies the line of the z_t_matrix by the theta_i_tH scalar

denominator = theta_1_tH' * z_t(IRF_hor:end)';

IRF = (1/denominator) * numerator;
IRF = IRF';

end


