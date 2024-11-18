clc; clear; close all;

%%%%%
% Script to use Monte Carlo simulation to test the validity of the dependent 
% wild bootstrapping for the LP and SVAR model average IRF estimator

% Use a VAR(1) model as the DGP
%%%%%

% Set random seed for reproducibility
rng(0);

% Simulation parameters
n_MC = 1000;          % Number of Monte Carlo simulations
T_total = 400;
T_burn = 100;
n_straps = 400;
clevel = 1.96;        % we are interested 95% CI
H = 20;               % IRF horizon
n_vars = 3;           % Number of variables
rind = 2;             % response variable index
sind = 1;             % shock variable index
coverage = zeros(n_MC, H);
nlag = 2;
bootstrap_ci_method = "normal";
use_bootstrap = 0;

% VAR(p) DGP
[B, AA, ~, optlag] = generate_DGP_rz(10,0);

%% Get true IRFs

for i = 1:length(AA)
    AA{i} = tril(AA{i});
end

for i = 1:length(AA)
    eval(['A' num2str(i) '= AA{' num2str(i) '};']);
end

% Structural shock variances (diagonal matrix)
D = eye(3); 

% Construct companion matrix for VAR(4)
A_toprow = AA{1};
for nlag=2:optlag
    A_toprow = [A_toprow, AA{nlag}];
end

A_rest = [];
for nlag=1:optlag-1
    A_rest_row = [zeros(n_vars, n_vars*(nlag-1)), eye(n_vars), zeros(n_vars, n_vars*(optlag-nlag))];
    A_rest = [A_rest; A_rest_row];
end

A = [A_toprow; A_rest];
Binv = inv(B);

% Compute true structural IRFs
true_irf = zeros(H, 1);
impact = Binv(:,sind);  % Initial impact of structural shock sind
true_irf(1) = impact(rind);  % Store response of variable rind

disp(eig(A))
disp(['Max eigenvalue is ' num2str(max(abs(eig(A))))])

% Compute subsequent horizons for true IRF
state = [impact; zeros((optlag-1)*n_vars,1)]; % State vector includes 4 lags of 3 variables
for h = 2:H
    state = A * state;
    true_irf(h) = state(rind);
end

figure;
plot(true_irf)
title('True IRF')

%% LP Experiment - Find Bias and Variance and Save

bias_store = zeros(n_MC, H);
irf_store = zeros(n_MC,H);
irf_squared_dev_true = zeros(n_MC,H);

for i = 1:n_MC
    disp(['Running MC rep number ', num2str(i)]);
    Y_total = zeros(T_total, 3);
    eps = mvnrnd(zeros(3,1), D, T_total);
    for t = optlag+1:T_total
        By_t = eps(t,:);
        for lag = 1:optlag
            By_t = By_t + Y_total(t-lag,:) * eval(['A' num2str(lag)]);
        end
        Y_total(t,:) = By_t * inv(B);
    end
    
    Y = Y_total(T_burn+1:end,:);

    % LP
    yresp = Y(:,rind);
    shock = Y(:,sind);
    ylag = lagmatrix(Y, [1:nlag]);
    ylag = ylag(nlag+1:end,:);
    Teff = size(ylag, 1);
    
    % create final Y and X data vectors for LP
    x = [ones(Teff,1) shock(nlag+1:end,:) ylag];
    yresp = yresp(nlag+1:end,:);

    rpos = 2;
    transformation = 1;
    opt = 0;

    % Estimate linlp on MC sample
    [liny, ~] = linlp(yresp, x, H, rpos, transformation, clevel, opt, use_bootstrap, nlag, n_straps, bootstrap_ci_method);
    
    irf_store(i,:) = liny;
    bias_store(i,:) = liny - true_irf';
    irf_squared_dev_true(i,:) = (liny - true_irf').^2;
end

% Results
var_irf = var(irf_store,1,1);
% var_irf = mean(irf_squared_dev_true,1);
std_irf = var_irf.^(1/2);
abs_meanbias = abs(mean(bias_store,1));

res = struct();
res.var = var_irf;
res.std = std_irf;
res.bias = abs_meanbias;
res.rmse = mean(irf_squared_dev_true,1).^(1/2);

filename = "Bootstrap/results_biasvar/rzdgp_lp_MC" + num2str(n_MC) + "_BS" + num2str(n_straps) + ".mat";
save(filename, "res")

%% SVAR Experiment - Find Bias and Variance and Save

bias_store = zeros(n_MC, H);
irf_store = zeros(n_MC,H);
irf_squared_dev_true = zeros(n_MC,H);

for i = 1:n_MC
    disp(['Running MC rep number ', num2str(i)]);
    eps = mvnrnd(zeros(3,1), D, T_total);

    Y_total = zeros(T_total, 3);
    eps = mvnrnd(zeros(3,1), D, T_total);
    for t = optlag+1:T_total
        By_t = eps(t,:);
        for lag = 1:optlag
            By_t = By_t + Y_total(t-lag,:) * eval(['A' num2str(lag)]);
        end
        Y_total(t,:) = By_t * inv(B);
    end
    
    Y = Y_total(T_burn+1:end,:);

    % VAR specs
    normalise = 1;
    y = Y;
    
    % Estimate SVAR on MC sample
    [liny, ~] = linSVAR(y,H,nlag,rind,sind,clevel,use_bootstrap,n_straps,normalise,bootstrap_ci_method);

    irf_store(i,:) = liny;
    bias_store(i,:) = liny - true_irf';
    irf_squared_dev_true(i,:) = (liny - true_irf').^2;

end

% Results
var_irf = var(irf_store,1,1);
% var_irf = mean(irf_squared_dev_true,1);
std_irf = var_irf.^(1/2);
abs_meanbias = abs(mean(bias_store,1));

res = struct();
res.var = var_irf;
res.std = std_irf;
res.bias = abs_meanbias;
res.rmse = mean(irf_squared_dev_true,1).^(1/2);

filename = "Bootstrap/results_biasvar/rzdgp_SVAR_MC" + num2str(n_MC) + "_BS" + num2str(n_straps) + ".mat";
save(filename, "res")

%% Mdl Avg Experiment 0.5 - Find Bias and Variance and Save

lambda = 0.5 * ones(1,H)

bias_store = zeros(n_MC, H);
irf_store = zeros(n_MC,H);
irf_squared_dev_true = zeros(n_MC,H);

for i = 1:n_MC
    iterStart = tic;
    disp(['Running MC rep number ', num2str(i)]);
    % Generate data from VAR(1)
    Y_total = zeros(T_total, 3);
    eps = mvnrnd(zeros(3,1), D, T_total);
    for t = optlag+1:T_total
        By_t = eps(t,:);
        for lag = 1:optlag
            By_t = By_t + Y_total(t-lag,:) * eval(['A' num2str(lag)]);
        end
        Y_total(t,:) = By_t * Binv;
    end
    
    Y = Y_total(T_burn+1:end,:);

    % LP
    yresp = Y(:,rind);
    shock = Y(:,sind);
    ylag = lagmatrix(Y, [1:nlag]);
    ylag = ylag(nlag+1:end,:);
    Teff = size(ylag, 1);
    
    % create final Y and X data vectors for LP
    x = [ones(Teff,1) shock(nlag+1:end,:) ylag];
    yresp = yresp(nlag+1:end,:);

    rpos = 2;
    transformation = 1;
    opt = 0;

    % VAR specs
    normalise = 1;
    y = Y;
    
    % Estimate linSVARLP_avg on MC sample
    [liny, ~] = linSVARLP_avg(H,nlag,bootstrap_ci_method, clevel, n_straps, lambda, use_bootstrap, ...
        yresp,x,rpos,transformation,opt, ...
        y,rind,sind,normalise);

    irf_store(i,:) = liny;
    bias_store(i,:) = liny - true_irf';
    irf_squared_dev_true(i,:) = (liny - true_irf').^2;
end

% Results
var_irf = var(irf_store,1,1);
% var_irf = mean(irf_squared_dev_true,1);
std_irf = var_irf.^(1/2);
abs_meanbias = abs(mean(bias_store,1));

res = struct();
res.var = var_irf;
res.std = std_irf;
res.bias = abs_meanbias;
res.rmse = mean(irf_squared_dev_true,1).^(1/2);

filename = "Bootstrap/results_biasvar/rzdgp_SVARLPavg0.5_MC" + num2str(n_MC) + "_BS" + num2str(n_straps) + ".mat";
save(filename, "res")

%% Mdl Avg Experiment Opt - Find Bias and Variance and Save

lambda_read = load("Weight Selection/opt_weight_v2.mat");
lambda = lambda_read.opt_weights

bias_store = zeros(n_MC, H);
irf_store = zeros(n_MC,H);
irf_squared_dev_true = zeros(n_MC,H);

for i = 1:n_MC
    disp(['Running MC rep number ', num2str(i)]);
    % Generate data from VAR(1)
    Y_total = zeros(T_total, 3);
    eps = mvnrnd(zeros(3,1), D, T_total);
    for t = optlag+1:T_total
        By_t = eps(t,:);
        for lag = 1:optlag
            By_t = By_t + Y_total(t-lag,:) * eval(['A' num2str(lag)]);
        end
        Y_total(t,:) = By_t * Binv;
    end
    
    Y = Y_total(T_burn+1:end,:);

    % LP
    yresp = Y(:,rind);
    shock = Y(:,sind);
    ylag = lagmatrix(Y, [1:nlag]);
    ylag = ylag(nlag+1:end,:);
    Teff = size(ylag, 1);
    
    % create final Y and X data vectors for LP
    x = [ones(Teff,1) shock(nlag+1:end,:) ylag];
    yresp = yresp(nlag+1:end,:);

    rpos = 2;
    transformation = 1;
    opt = 0;

    % VAR specs
    normalise = 1;
    y = Y;
    
    % Estimate linSVARLP_avg on MC sample
    [liny, ~] = linSVARLP_avg(H,nlag,bootstrap_ci_method, clevel, n_straps, lambda, use_bootstrap, ...
        yresp,x,rpos,transformation,opt, ...
        y,rind,sind,normalise);

    irf_store(i,:) = liny;
    bias_store(i,:) = liny - true_irf';
    irf_squared_dev_true(i,:) = (liny - true_irf').^2;
end

% Results
var_irf = var(irf_store,1,1);
% var_irf = mean(irf_squared_dev_true,1);
std_irf = var_irf.^(1/2);
abs_meanbias = abs(mean(bias_store,1));

res = struct();
res.var = var_irf;
res.std = std_irf;
res.bias = abs_meanbias;
res.rmse = mean(irf_squared_dev_true,1).^(1/2);

filename = "Bootstrap/results_biasvar/rzdgp_SVARLPavgoptv2_MC" + num2str(n_MC) + "_BS" + num2str(n_straps) + ".mat";
save(filename, "res")

%% Plot

plot_files = {
    "Bootstrap/results_biasvar/rzdgp_SVARLPavgoptv2_MC1000_BS400.mat",
    "Bootstrap/results_biasvar/rzdgp_SVARLPavg0.5_MC1000_BS400.mat",
    "Bootstrap/results_biasvar/rzdgp_SVAR_MC1000_BS400.mat",
    "Bootstrap/results_biasvar/rzdgp_lp_MC1000_BS400.mat"
};

% Load data
for i=1:length(plot_files)
    eval(['res' num2str(i) '= load(plot_files{i});'])
    eval(['res' num2str(i) '= res' num2str(i) '.res;'])
    eval(['bias' num2str(i) '=res' num2str(i) '.bias;'])
    eval(['std' num2str(i) '=res' num2str(i) '.std;'])
end

colors = {
    [0.8500 0.3250 0.0980], % Deep orange (avg)
    [0 0.4470 0.7410],      % Rich blue (avg)
    [0.3 0.3 0.3],          % Light gray (base)
    [0.6 0.6 0.6]           % Darker gray (base)
};

colors = ['g', 'b', 'k', 'm'];
% Create figure with reasonable size
fig = figure('Color', 'white', 'Position', [100 100 1000 400]); % Half the previous size

% Main title with improved formatting
sgtitle({'\bf Bias and Standard Deviation of Estimators', ...
    ['\rm MC Simulations: ', num2str(n_MC)]}, ...
    'FontSize', 11)

% Bias plot
subplot(1,2,1)
title('\bf Average Absolute Bias', 'FontSize', 10)
hold on
for i=1:length(plot_files)
    eval(['bias_data = bias' num2str(i) ';'])
    plot(0:H-1, bias_data, 'LineStyle', '-', 'LineWidth', 2, 'Color', colors(i));
end
legend('Model Avg. Opt', 'Model Avg. 0.5', 'SVAR', 'LP', ...
    'Location', 'northeast', 'FontSize', 10, 'Box', 'on')
xlabel('Horizon', 'FontSize', 9)
ylabel('Bias', 'FontSize', 9)
grid on
box on
hold off

% Variance plot
subplot(1,2,2)
title('\bf Standard Deviation', 'FontSize', 10)
hold on
for i=1:length(plot_files)
    eval(['std_data = std' num2str(i) ';'])
    plot(0:H-1, std_data, 'LineStyle', '-', 'LineWidth', 2, 'Color', colors(i));
end
legend('Model Avg. Opt', 'Model Avg. 0.5', 'SVAR', 'LP', ...
    'Location', 'southeast', 'FontSize', 10, 'Box', 'on')
xlabel('Horizon', 'FontSize', 9)
ylabel('Standard Deviation', 'FontSize', 9)
grid on
box on
hold off

% Add spacing between subplots
fig.Position = [0,0,800, 300]

% Save the figure
saveas(fig, 'Fig/bias_std_rzdgp_allest.png')

%% RMSE

% Load data
for i=1:length(plot_files)
    eval(['res' num2str(i) '= load(plot_files{i});'])
    eval(['res' num2str(i) '= res' num2str(i) '.res;'])
    eval(['rmse' num2str(i) '=res' num2str(i) '.rmse;'])
end

colors = {
    [0.8500 0.3250 0.0980], % Deep orange (avg)
    [0 0.4470 0.7410],      % Rich blue (avg)
    [0.3 0.3 0.3],          % Light gray (base)
    [0.6 0.6 0.6]           % Darker gray (base)
};

colors = ['g', 'b', 'k', 'm'];


% Create figure with reasonable size
fig = figure('Color', 'white', 'Position', [100 100 1000 400]); % Half the previous size


title(['RMSE of Estimators',' | MC Simulations: ', num2str(n_MC)])
hold on
for i=1:length(plot_files)
    eval(['rmse_data = rmse' num2str(i) ';'])
    plot(0:H-1, rmse_data, 'LineStyle', '-', 'LineWidth', 2, 'Color', colors(i));
end
legend('Model Avg. Opt', 'Model Avg. 0.5', 'SVAR', 'LP', ...
    'Location', 'northeast', 'FontSize', 10, 'Box', 'on')
xlabel('Horizon', 'FontSize', 9)
ylabel('Bias', 'FontSize', 9)
grid on
box on
hold off

% Add spacing between subplots
set(gcf, 'Units', 'normalized');
p = get(fig, 'Position');
p(3) = p(3) * 1.05; % Just slight width increase for spacing
set(fig, 'Position', p);

% Save the figure
saveas(fig, 'Fig/rmse_rzdgp_allest.png')
