function [IRF,n_lags_est] = LP_GLS_est(data_sim,settings);
% Function for estimating IRFs using generalised least squares LP using
% Lusompa method of transforming the dependent variable

% preparations

run('Estimation_Setup'); % common setup for all estimation methods

% estimate IRF via LP
[Beff, Sigma_hold] = lp_gls_analytical(Y', nlags, IRF_hor, 0);

IRF = arrayfun(@(i) Beff(end, responseV, i, 5), 1:IRF_hor);

% normalize - TBD
% IRF_normalize = IRF_LP(Y,recursiveShock,normalizeV,nlags,0);

%
IRF = IRF';

end

