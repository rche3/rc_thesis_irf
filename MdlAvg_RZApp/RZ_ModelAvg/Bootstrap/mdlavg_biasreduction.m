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

%% Mdl Avg Experiment 0.5 - Find Bias and Variance and Save

weights = [0:0.05:1];
abs_bias_per_weight = zeros(size(weights,2),H);

for w = 1:size(weights,2)
    weight = weights(w);
    lambda = weight * ones(1,H);
    
    bias_store = zeros(n_MC, H);
    
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
    
        bias_store(i,:) = liny - true_irf';
    end
    abs_meanbias = abs(mean(bias_store,1));
    abs_bias_per_weight(w,:) = abs_meanbias;
end

%% Calculate best weights

sumhorizon_abs_bias_per_weight = sum(abs_bias_per_weight, 2)

[sorted_values, indices] = sort(sumhorizon_abs_bias_per_weight);
lowest_indices = indices(1:5)

%% Plot

plot_indices = [lowest_indices', 21]; % plot lowest 3 and SVAR (w=1)
fig = figure;

hold on
for w=1:length(plot_indices)
    plot_idx = plot_indices(w);
    data_plot = abs_bias_per_weight(plot_idx,:);
    plot(1:H, data_plot, 'LineWidth', 1.2)
end

legendinput = cell(1,length(plot_indices));
for i = 1:length(plot_indices)
    legendinput{i} = num2str(weights(plot_indices(i)));
end
legend(legendinput)
xlabel('Horizon')
ylabel('Absolute Mean Bias')

saveas(fig, 'Fig/bias_reduction_mdlavg.png')