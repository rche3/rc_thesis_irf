
clear all;

%% init conditions
cd /Users/rogerchen/Documents/MATLAB/rc_thesis_irf/DFM/

%%%% spec variables
spec_id = 1; % seed for random draws of specifications (= DGPs from encompassing model)
lag_type = 4; % No. of lags to impose in estimation, or NaN (= AIC)
mode_type = 1; % robustness check mode:
               % 1 (baseline), 2 (small sample), 3 (large sample),
               % 4 (salient series), 5 (more observables), 6 (first diff)
estim_diagn = 0; % =1: show DFM estimation diagnostics

%%% SPEC TYPES
dgp_type = 'G'; % structural shock: either 'G' or 'MP'
estimand_type = 'ObsShock'; % structural estimand: either 'ObsShock', 'Recursive', or 'IV'
data = load('sample_data_obsG.mat');
sample_data = data.data_sim_select;

%%% SETTINGS
run(fullfile('Settings', 'shared'));
run(fullfile('Settings', dgp_type));
run(fullfile('Settings', estimand_type));
run(fullfile('Settings', 'check_mode'));

[irf1, ~] = LP_est(sample_data, settings, 0);
[irf2, ~] = SVAR_est(sample_data, settings, 0);
[irf3, ~] = LP_GLS_est(sample_data, settings);

horizon_0_response_OLS = irf1(1)
horizon_0_response_GLS = irf3(1)

%%% plot functions
figure;
ax1 = subplot(1,4,1);
plot(irf1)
title('IRF with OLS LP')

ax2 = subplot(1,4,2);
plot(irf2)
title('IRF with SVAR')

ax3 = subplot(1,4,3);
plot(irf3)
title('IRF with GLS LP')

% linkaxes([ax1, ax2], 'y')
linkaxes([ax1, ax2, ax3], 'y')
% linkaxes([ax1, ax2, ax3, ax4], 'y')
