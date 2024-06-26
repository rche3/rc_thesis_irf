function [IRF,n_lags_est] = LP_GLS_est(data_sim,settings);
% Function for estimating IRFs using least-squares LP

% preparations

run('Estimation_Setup'); % common setup for all estimation methods

% estimate IRF via LP
blocksize = 20;
bootstrap_reps = 1000;

[Beff, Sigma_hold] = lusompa_lp_gls_boot(Y', nlags, bootstrap_reps, IRF_hor, blocksize, 0)

IRF = zeros(IRF_hor, 1);

for i=1:IRF_hor
    IRF(i) = Beff(end, responseV, i, 1000)
end

% normalize - TBD
% IRF_normalize = IRF_LP(Y,recursiveShock,normalizeV,nlags,0);

%
IRF = IRF';

end

