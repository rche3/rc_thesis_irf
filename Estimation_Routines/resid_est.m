function [IRF, n_lags_est] = resid_est(data_sim,settings)
%RESID_EST Summary of this function goes here
%   Detailed explanation goes here

run('Estimation_Setup'); % common setup for all estimation methods

%%% compute the LP residuals (LP errors)
res = IRF_LP_res(Y,recursiveShock,responseV,nlags,IRF_hor - 1); % residuals for IRF estimator
disp(size(res))
disp(res)

% placeholders
IRF = 1
n_lags_est = 1
end


