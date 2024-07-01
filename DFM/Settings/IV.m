%% SPECIFIC SETTINGS FOR IV IDENTIFICATION

% IV equation

DF_model.IV.rho     = 0.25; % baseline IV persistence
DF_model.IV.alpha   = 1; % IV shock coefficient
DF_model.IV.sigma_v = 1; % baseline IV noise

settings.est.IV.IV_persistence_scale = [0 1 2]; % scale up or down baseline persistence by this factor (use 1 if only want baseline persistence)
settings.est.IV.IV_strength_scale = [1.1 1.5 2.3]; % scale up or down baseline IV strength by this factor (use 1 if only want baseline strength)

% estimation methods

settings.est.methods_name    = [settings.est.methods_name {'svar_iv'} {'resid_est'} {'lp_IV_controls'}]; % add SVAR-IV, resid_est, ... (other IV estimators) to the estimator list

% indicate estimand

settings.est.with_shock      = 0; % shock is observed and ordered first in data?
settings.est.recursive_shock = 0; % use recursive shock
settings.est.with_IV         = 1; % IV is ordered first in data or used in IV method?