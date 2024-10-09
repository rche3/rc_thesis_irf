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
Beff_normalise = lp_gls_analytical_simp(Y', nlags, 1, recursiveShock, recursiveShock);
IRF_normalize_gls = Beff_normalise(1)
% IRF = IRF/IRF_normalize_gls;

IRF = IRF';

end

