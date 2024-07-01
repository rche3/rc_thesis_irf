function [IRF,n_lags_est] = LP_IV_controls_est(data_sim,settings)
%LP_IV_CONTROLS_EST Summary of this function goes here

% preparations

run('Estimation_Setup'); % common setup for all estimation methods

% note recursiveShock = 1 for this estimator since this is the definition
% of IV and observed shock.

% Y is a T x 6 vector with the IV in the first column

[IRF,w] = IRF_LP_IVC(Y,recursiveShock,responseV,nlags,IRF_hor - 1); % IRF to one unit of shock

end

