function [IRF,n_lags_est] = LP_GLS_RC_est(data_sim,settings);
% Function for estimating IRFs using GLS LP using
% Lusompa method of transforming the dependent variable

% preparations
run('Estimation_Setup'); % common setup for all estimation methods

% estimate IRF via GLS LP function, first part is regressing y_{t+h} on z_t
fgls_method = 2; % 1 = Lusompa, 2 = Breitung
Beff = lp_gls_analytical_manual(Y, nlags, IRF_hor, recursiveShock, responseV, fgls_method);
IRF = Beff;

% need to divide this by the regression of i_t on z_t
IRF_normalize_gls = lp_gls_analytical_manual(Y, nlags, 0, recursiveShock, normalizeV, fgls_method);

IRF = IRF/IRF_normalize_gls;

% disp(['GLS result is ', num2str(IRF(1))]);

end

