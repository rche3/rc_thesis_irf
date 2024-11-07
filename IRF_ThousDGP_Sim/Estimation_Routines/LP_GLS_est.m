function [IRF,n_lags_est] = LP_GLS_est(data_sim,settings);
% Function for estimating IRFs using GLS LP using
% Lusompa method of transforming the dependent variable

% preparations
run('Estimation_Setup'); % common setup for all estimation methods

% estimate IRF via LP
Beff = lp_gls_analytical_simp(Y', nlags, IRF_hor, recursiveShock, responseV);

% IRF = squeeze(Beff(responseV, responseV, 1:IRF_hor));
IRF = Beff;

% normalize?
IRF_normalize_gls = lp_gls_analytical_simp(Y', nlags, 0, recursiveShock, normalizeV);
IRF = IRF/IRF_normalize_gls;

% disp(['GLS result is ', num2str(IRF(1))]);

end

