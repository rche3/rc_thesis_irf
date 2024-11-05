clc
clear
close all

rng(123456)

sim_name = 'ReviewPaperResults_all_2021_05_06_';

%% Set Simulation Parameters

nsim            = 1000;
nboot           = 2000;
T_vec           = [100;500];    % 200,500
a11_vec         = [1;0.95;0.9;0.5;0.1];   % [1;0.95;0.9;0.5;0.1];
sigma2_nu_vec   = [0.2346;3];      % [0.2346];
p_vec           = [1;12];      % 12
d_vec           = [1];
H_vec           = [20];     % 20
bias_correction = [1];        % 0

total_nb_spec = length(T_vec)*length(a11_vec)*length(d_vec)*length(sigma2_nu_vec)*length(p_vec)*length(H_vec)

p_true      = 1;
B           = [1 0; 0.5 3];
k           = 2;
psi         = 1;

sigma2_w    = eye(2);
% y0          = zeros(2,1)
alpha_boot  = 10; % bootstrap quantiles (set to 10 for 90% CI)
boot_method = 6;    % 1 residual resample
                    % 2 wild with normal
                    % 3 wild with Rademacher
                    % 4 uses new draw of U_true, z_true
                    % 5 test run, similar to SW2017
                    % 6 Moving Block Bootstrap
                    % 7 new bootstrap
                    
boot_names = {'','Wild','','','','MBB',''};                    

length_type     = 'average';
% length_type     = 'median';

curr_dir = pwd;
dir_case            = [sim_name,num2str(boot_names{boot_method}), '_bc',num2str(bias_correction)]; % where to save the simulation results (not the figures) 
mkdir(dir_case)      
cd(dir_case)
mkdir('WorkSpace')
mkdir('Figures')
cd(curr_dir)   

const       = 1;    % 1: The VAR estimates a constant (although there is none in the DGP)
                    % 0: No constant in VAR, but data is demeaned 

const_proxy_eq      = 1;    % 1: The Poxy Equation includes a constant
                            % 0: No constant in Proxy Equation                   
                    
count_sim       = 0;
time_sim_all    = NaN(total_nb_spec,1);
case_name_all   = cell(total_nb_spec,1);

for TT = 1:length(T_vec)
    T       = T_vec(TT);
    for aa = 1:length(a11_vec)
        a11     = a11_vec(aa);
        A1      = [a11 0 ; 0.5 0.5];
        for dd = 1:length(d_vec)
            d   = d_vec(dd);
            for ss = 1:length(sigma2_nu_vec)
                sigma2_nu = sigma2_nu_vec(ss);
                for pp = 1:length(p_vec)
                    p   = p_vec(pp);
                    for HH = 1:length(H_vec)
                        H = H_vec(HH);
% T   = 100;
% a11 = 0.9;
% A1  = [a11 0 ; 0.5 0.5];
% d   = 1;
% sigma2_nu = 0.2346;
% p   = 1;
% H = 20;

corr_wz = psi*sqrt(d)*sqrt(sigma2_w(1,1))/sqrt(psi^2*sigma2_w(1,1)+sigma2_nu);

IRF_true    = zFC_IRFs(A1, B./B(1,1), H+1, 1, 1);

case_name = ['N', num2str(nsim),...
        '_B', num2str(nboot),...
        '_T', num2str(T),...
        '_a', num2str(a11),...
        '_d', num2str(d),...
        '_CorrWZ', num2str(round(corr_wz,1)),...
        '_p', num2str(p),...
        '_H', num2str(H),...
        '_BC', num2str(bias_correction),...
        '_Boot', boot_names{boot_method}];
    
[coverage_var, av_length_var,mse_var,...
          coverage_lp, av_length_lp,mse_lp,...
          coverage_lp_aug, av_length_lp_aug,mse_lp_aug,...
          coverage_lp_iv, av_length_lp_iv,mse_lp_iv,...
          coverage_lp_iv_c_no_z, av_length_lp_iv_c_no_z,mse_lp_iv_c_no_z,...
          coverage_lp_iv_c, av_length_lp_iv_c,mse_lp_iv_c,...
          coverage_new, av_length_new,mse_new,...
          coverage_new_ss, av_length_new_ss,mse_new_ss,...
          coverage_bb, av_length_bb,mse_bb,...
          coverage_bb_gls, av_length_bb_gls,mse_bb_gls,...
          coverage_lus_gls, av_length_lus_gls,mse_lus_gls,...
          time_sim]        = zFC_SimulateSurvey_04(A1,B,d ,sigma2_w,...
                                                            psi, sigma2_nu, p_true,T,nsim,...
                                                                H,p,const,const_proxy_eq,k,...
                                                                nboot, alpha_boot,...
                                                                IRF_true, boot_method,bias_correction,length_type);           
     

curr_dir = pwd;
cd(dir_case)
cd('WorkSpace')
save(['WS_', case_name, '.mat'], 'nsim','nboot', 'T', 'a11', 'H', 'd', 'corr_wz', 'p',...
         'coverage_var', 'av_length_var','mse_var',...
          'coverage_lp', 'av_length_lp','mse_lp',...
          'coverage_lp_aug', 'av_length_lp_aug','mse_lp_aug',...
          'coverage_lp_iv', 'av_length_lp_iv','mse_lp_iv',...
          'coverage_lp_iv_c_no_z', 'av_length_lp_iv_c_no_z','mse_lp_iv_c_no_z',...
          'coverage_lp_iv_c', 'av_length_lp_iv_c','mse_lp_iv_c',...
          'coverage_new', 'av_length_new','mse_new',...
          'coverage_new_ss', 'av_length_new_ss','mse_new_ss',...
          'coverage_bb', 'av_length_bb','mse_bb',...
          'coverage_bb_gls', 'av_length_bb_gls','mse_bb_gls',...
          'coverage_lus_gls', 'av_length_lus_gls','mse_lus_gls',...
            'time_sim',...
    'boot_method','bias_correction','length_type',...
    'p_true','IRF_true', 'B', 'k', 'psi', 'sigma2_w', 'sigma2_nu')

cd(curr_dir)                          

count_sim = count_sim + 1;
time_sim_all(count_sim) = time_sim;
case_name_all{count_sim} = case_name;
disp([num2str(total_nb_spec-count_sim)])                 

                    end
                end
            end
        end
    end
end

total_running_time_in_minutes = sum(time_sim_all)/60;

disp('DONE')

