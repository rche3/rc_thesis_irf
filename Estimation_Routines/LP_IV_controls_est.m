function [IRF,n_lags_est] = LP_IV_controls_est(data_sim,settings)
% preparations
run('Estimation_Setup'); % common setup for all estimation methods

% note recursiveShock = 1 for this estimator since this is the definition of IV and observed shock.
% Y is a T x 6 vector with the IV in the first column
% Y1 corresponds to the policy variable <-> i_t (Gov spend or Monetary policy)

IRF = zeros(IRF_hor, 1);

% loop over horizons:
for h = 0:IRF_hor-1
    % Step 1: Obtain the "projection residuals" Y_t^\perp and z_t^\perp (Stock and Watson, 2018)
    yr_res = proj_w_resid(responseV, Y, h, recursiveShock, nlags, h, include_proxy);
    z_res = proj_w_resid(recursiveShock, Y, 0, recursiveShock, nlags, h, include_proxy);
    y1_res = proj_w_resid(policyV, Y, 0, recursiveShock, nlags, h, include_proxy); % this is the location of the policy variable
    
    % Step 2: Compute the estimator as the fraction of the two sums
    num = sum(yr_res.*z_res);
    denom = sum(y1_res.*z_res);
    IRF(h+1) = num / denom;
end

end

