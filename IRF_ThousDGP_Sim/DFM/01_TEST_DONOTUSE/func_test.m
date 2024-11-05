%% USED TO TEST G AND IV estimation
clear all;,

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
estimand_type = 'IV'; % structural estimand: either 'ObsShock', 'Recursive', or 'IV'
data = load('sample_data_IVg.mat');
sample_data = data.sample_data_IVg;

run(fullfile('Settings', 'shared'));
run(fullfile('Settings', dgp_type));
run(fullfile('Settings', estimand_type));
run(fullfile('Settings', 'check_mode'));

[irf1, ~] = LP_est(sample_data, settings, 0);
[irf2, ~] = LP_IV_controls_est(sample_data, settings);
[irf3, ~] = LP_GLS_est(sample_data, settings);
[irf4, ~] = LP_GLS_RC_est(sample_data, settings);

%%
[irf5, ~] = resid_est(sample_data, settings);

% plot functions
fig = figure;
ax1 = subplot(1,5,1);
plot(irf1)
title('IRF with OLS LP')
axis tight

ax2 = subplot(1,5,2);
plot(irf2)
title('IRF with LP-IV')
axis tight

ax3 = subplot(1,5,3);
plot(irf3)
title('IRF with GLS LP')
axis tight

ax4 = subplot(1,5,4);
plot(irf4)
title('IRF with manual GLS LP code')
axis tight

ax5 = subplot(1,5,5);
plot(irf5)
title('IRF with resid est')
axis tight

% linkaxes([ax1, ax2], 'y')
% linkaxes([ax1, ax2, ax4], 'y')
linkaxes([ax1, ax2, ax3, ax4 ax5], 'y');

fig.Position = [100, 100, 1400, 400];
