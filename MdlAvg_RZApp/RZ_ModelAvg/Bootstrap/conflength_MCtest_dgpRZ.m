clc; clear; close all;

%%%%%
% Script to use Monte Carlo simulation to test the validity of the dependent 
% wild bootstrapping for the LP and SVAR model average IRF estimator

% Use a VAR(1) model as the DGP
%%%%%

% Set random seed for reproducibility
rng(42);

% Simulation parameters
n_MC = 100;          % Number of Monte Carlo simulations
T_total = 400;
T_burn = 100;
H = 20;               % IRF horizon
n_vars = 3;           % Number of variables
rind = 2;             % response variable index
sind = 1;             % shock variable index
coverage = zeros(n_MC, H);

% VAR(p) DGP
[B, AA, ~, optlag] = generate_DGP_rz(10,0);

%%
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

%% GENERAL SETTINGS
n_straps = 400;
n_MC = 1000;

%% Monte Carlo simulation - LP
% specs for model average estimation
nlag = 2;
bootstrap_ci_method = "normal";
use_bootstrap = 1;
clevel = 1.96; % we are interested 95% CI
timePerLoop = zeros(n_MC, 1);  % Store time for each iteration


ci_length = zeros(n_MC, H);

for i = 1:n_MC
    % Generate data from VAR(1)
    iterStart = tic;  % Start timing this iteration

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
    [liny, confidencey, ~, ~, ~] = linlp(yresp, x, H, rpos, transformation, clevel, opt, use_bootstrap, nlag, n_straps, bootstrap_ci_method);
    
    lower_ci = confidencey(1,:,:);
    upper_ci = confidencey(2,:,:);

    ci_length(i,:) = upper_ci - lower_ci;
    
    % Calculate and display timing for this iteration
    timePerLoop(i) = toc(iterStart);
    estimatedTimeRemaining = mean(timePerLoop(1:i)) * (n_MC - i);
    
    fprintf('Iteration %d/%d completed in %.2f seconds\n', i, n_MC, timePerLoop(i));
    fprintf('Estimated time remaining: %.2f minutes\n', estimatedTimeRemaining/60);
end

% save
filename = "Bootstrap/results_conflength/rzDGP_lp_conflength_MC" + num2str(n_MC) + "_BS" + num2str(n_straps) + ".mat";
save(filename, "ci_length")

%% Monte Carlo simulation - SVAR with bootstrap

% specs for model average estimation

clear ci_length
nlag = 2;
use_bootstrap = 1;
bootstrap_ci_method = "normal";
clevel = 1.96; % we are interested 95% CI
timePerLoop = zeros(n_MC, 1);  % Store time for each iteration

ci_length = zeros(n_MC, H);

for i = 1:n_MC
    % Generate data from VAR(1)
    iterStart = tic;  % Start timing this iteration
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
    [liny, confidencey] = linSVAR(y,H,nlag,rind,sind,clevel,use_bootstrap,n_straps,normalise,bootstrap_ci_method);
     
    ci_length(i,:) = upper_ci - lower_ci;

    % Timing
    timePerLoop(i) = toc(iterStart);
    estimatedTimeRemaining = mean(timePerLoop(1:i)) * (n_MC - i);
    
    fprintf('Iteration %d/%d completed in %.2f seconds\n', i, n_MC, timePerLoop(i));
    fprintf('Estimated time remaining: %.2f minutes\n', estimatedTimeRemaining/60);
end

% save
filename = "Bootstrap/results_conflength/rzDGP_svar_conflength_MC" + num2str(n_MC) + "_BS" + num2str(n_straps) + ".mat";
save(filename, "ci_length")

%% Monte Carlo simulation - model average
% settings
clear ci_length
nlag = 2;
bootstrap_ci_method = "normal";
clevel = 1.96; % we are interested 95% CI
% lambda_read = load("Weight Selection/opt_weight.mat");
% lambda = lambda_read.opt_weights
lambda = 0.5 * ones(1,H );
timePerLoop = zeros(n_MC, 1);  % Store time for each iteration

ci_length = zeros(n_MC, H);

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
    [liny, confidencey] = linSVARLP_avg(H,nlag,bootstrap_ci_method, clevel, n_straps, lambda, use_bootstrap, ...
        yresp,x,rpos,transformation,opt, ...
        y,rind,sind,normalise);
    
    lower_ci = confidencey(1,:,:);
    upper_ci = confidencey(2,:,:);
    
    ci_length(i,:) = upper_ci - lower_ci;

    % Calculate and display timing for this iteration
    timePerLoop(i) = toc(iterStart);
    estimatedTimeRemaining = mean(timePerLoop(1:i)) * (n_MC - i);
    
    fprintf('Iteration %d/%d completed in %.2f seconds\n', i, n_MC, timePerLoop(i));
    fprintf('Estimated time remaining: %.2f minutes\n', estimatedTimeRemaining/60);
end

% save
filename = "Bootstrap/results_conflength/rzDGP_mdlavg0.5_conflength_MC" + num2str(n_MC) + "_BS" + num2str(n_straps) + ".mat";
save(filename, "ci_length")

%% Monte Carlo simulation - model average
clear ci_length
% settings
nlag = 2;
bootstrap_ci_method = "normal";
clevel = 1.96; % we are interested 95% CI
lambda_read = load("Weight Selection/opt_weight.mat");
lambda = lambda_read.opt_weights
timePerLoop = zeros(n_MC, 1);  % Store time for each iteration

ci_length = zeros(n_MC, H);

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
    [liny, confidencey] = linSVARLP_avg(H,nlag,bootstrap_ci_method, clevel, n_straps, lambda, use_bootstrap, ...
        yresp,x,rpos,transformation,opt, ...
        y,rind,sind,normalise);
    
    lower_ci = confidencey(1,:,:);
    upper_ci = confidencey(2,:,:);
    
    ci_length(i,:) = upper_ci - lower_ci;

    % Calculate and display timing for this iteration
    timePerLoop(i) = toc(iterStart);
    estimatedTimeRemaining = mean(timePerLoop(1:i)) * (n_MC - i);
    
    fprintf('Iteration %d/%d completed in %.2f seconds\n', i, n_MC, timePerLoop(i));
    fprintf('Estimated time remaining: %.2f minutes\n', estimatedTimeRemaining/60);
end

% save
filename = "Bootstrap/results_conflength/rzDGP_mdlavgopt_conflength_MC" + num2str(n_MC) + "_BS" + num2str(n_straps) + ".mat";
save(filename, "ci_length")

%% PLOT FOR PAPER - All Three
plot_files = {
    "Bootstrap/results_conflength/rzDGP_lp_conflength_500_BS20.mat",
    "Bootstrap/results_conflength/rzDGP_svar_conflength_500_BS20.mat",
    "Bootstrap/results_conflength/rzDGP_mdlavg0.5_conflength_MC500_BS20.mat",
    "Bootstrap/results_conflength/rzDGP_mdlavgopt_conflength_MC500_BS20.mat"
    };

for i=1:length(plot_files)
    eval(['res' num2str(i) '= load(plot_files{i});']);
    eval(['length_data' num2str(i) '=res' num2str(i) '.ci_length;']);
end

fig = figure;
hold on
colors = ['g', 'b', 'k', 'm', 'c', 'y'];  % Basic MATLAB colors

for i=1:length(plot_files)
    eval(['length_data = length_data' num2str(i)]); 
    mean_length = mean(length_data,1); % mean along columns / collapse rows
    plot(1:H, mean_length, 'LineStyle', '-', 'LineWidth', 1, 'Color', colors(i));
end
xlim([0,H])
xlabel('Horizon of IRF')
ylabel('Average Confidence Interval Length')
title('Monte Carlo Simulation of Estimator Confidence Interval Lengths')
legend('LP', 'SVAR', 'Model Avg. 0.5', 'Model Avg. Opt');

subtitle(['Number of MC sims: ', num2str(n_MC), '; Number of bootstrap reps: ', num2str(n_straps)]);

saveas(fig, 'Fig/estimator_ci_length.png')
