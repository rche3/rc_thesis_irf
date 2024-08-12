function [IRF,n_lags_est,largest_root,LM_stat,LM_pvalue,Hausman_stat,Hausman_pvalue,Granger_stat,Granger_pvalue] = VARLP_avg_est(data_sim,settings,bias_corrected);
    disp('estimating model averaged now!')
    %   Detailed explanation goes here
    run('Estimation_Setup'); % common setup for all estimation methods

    [IRF_var,~,largest_root,LM_stat,LM_pvalue,Hausman_stat,Hausman_pvalue,Granger_stat,Granger_pvalue] ...
    = SVAR_est(data_sim,settings,bias_corrected);
    
    [IRF_lp,n_lags_est] = LP_est(data_sim,settings,bias_corrected);

    IRF = IRF_var.*varlp_est_weights + IRF_lp.*(1-varlp_est_weights);
end
