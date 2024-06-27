function [IRF, n_lags_est] = resid_est(data_sim,settings)
%RESID_EST Summary of this function goes here
%   Detailed explanation goes here

run('Estimation_Setup'); % common setup for all estimation methods

%%% compute the LP residuals (LP errors)
res = IRF_LP_res(Y,recursiveShock,responseV,nlags,IRF_hor - 1); % residuals for IRF estimator

%%% compute the estimator using the proxy z_t
% restrict the z_t to p+1:T-H
z_t = Y(nlags+1:end,1)';

theta_tH = res;

window_size = IRF_hor;
num_windows = length(z_t) - window_size +1;

z_t_matrix = fliplr(hankel(z_t(1:num_windows), z_t(num_windows:end)));

% multiply each row of z_t_matrix by the corresponding scalar resid
z_t_matrix_scaled = bsxfun(@times, theta_tH, z_t_matrix);

numerator = sum(z_t_matrix_scaled, 1);

denominator = theta_tH' * z_t(IRF_hor:end)';

IRF = (1/denominator) * numerator

end


