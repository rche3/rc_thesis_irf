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

%% Monte Carlo simulation - LP
% specs for model average estimation
nlag = 2;
bootstrap_ci_method = "normal";
use_bootstrap = 1;
clevel = 1.96; % we are interested 95% CI
timePerLoop = zeros(n_MC, 1);  % Store time for each iteration
n_straps = 400;
n_MC = 1000;

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
    [liny, confidencey] = linlp(yresp, x, H, rpos, transformation, clevel, opt, use_bootstrap, nlag, n_straps, bootstrap_ci_method);
    
    lower_ci = confidencey(1,:,:);
    upper_ci = confidencey(2,:,:);
    
    % Check coverage
    for h = 1:H
        coverage(i,h) = true_irf(h) >= lower_ci(:,h) & true_irf(h) <= upper_ci(:,h);
    end
    temp_coverage = mean(coverage(1:i,:), 1);
    
    disp(['Current coverage is: ', num2str(temp_coverage)])
    % Calculate and display timing for this iteration
    timePerLoop(i) = toc(iterStart);
    estimatedTimeRemaining = mean(timePerLoop(1:i)) * (n_MC - i);
    
    fprintf('Iteration %d/%d completed in %.2f seconds\n', i, n_MC, timePerLoop(i));
    fprintf('Estimated time remaining: %.2f minutes\n', estimatedTimeRemaining/60);
end

% Results
mean_coverage = mean(coverage,1); % 1 = rows, to "collapse" the rows

%% Save LP
filename = "Bootstrap/results/v2_rzDGPVAR2_lp_MC" + num2str(n_MC) + "_BS" + num2str(n_straps) + ".mat";
save(filename, "mean_coverage")

%% Monte Carlo simulation - SVAR with bootstrap

% specs for model average estimation
nlag = 2;
use_bootstrap = 1;
bootstrap_ci_method = "normal";
clevel = 1.96; % we are interested 95% CI
timePerLoop = zeros(n_MC, 1);  % Store time for each iteration
n_straps = 400;
n_MC = 1000;

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

%     var1 = varm(3,nlag);
% 
%     estvar1 = estimate(var1,y);
%     
%     [response, lower, upper] = irf(estvar1);
%     
%     liny = response(:,sind,rind);
%     confidenceysvar(1,:,:) = lower(:,sind,rind);
%     confidenceysvar(2,:,:) = upper(:,sind,rind);
%     
    lower_ci = confidencey(1,:,:);
    upper_ci = confidencey(2,:,:);
    
    % Check coverage
    for h = 1:H
        coverage(i,h) = true_irf(h) >= lower_ci(:,h) & true_irf(h) <= upper_ci(:,h);
    end
    temp_coverage = mean(coverage(1:i,:), 1);
    
    disp(['Current coverage is: ', num2str(temp_coverage)])

    % Timing
    timePerLoop(i) = toc(iterStart);
    estimatedTimeRemaining = mean(timePerLoop(1:i)) * (n_MC - i);
    
    fprintf('Iteration %d/%d completed in %.2f seconds\n', i, n_MC, timePerLoop(i));
    fprintf('Estimated time remaining: %.2f minutes\n', estimatedTimeRemaining/60);
end

% Results
mean_coverage = mean(coverage,1); % 1 = rows, to "collapse" the rows

%% Save SVAR
filename = "Bootstrap/results/rzDGPVAR2_SVAR_MC" + num2str(n_MC) + "_BS" + num2str(n_straps) + ".mat";

save(filename, "mean_coverage")

%% Monte Carlo simulation - model average
% settings
nlag = 2;
bootstrap_ci_method = "normal";
clevel = 1.96; % we are interested 95% CI
lambda_read = load("Weight Selection/opt_weight.mat");
lambda = lambda_read.opt_weights
lambda = 0.5 * ones(1,H )
timePerLoop = zeros(n_MC, 1);  % Store time for each iteration
n_MC = 1000;
n_straps = 400;

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
    [liny, confidencey] = linSVARLP_avg(H,nlag,bootstrap_ci_method, clevel, n_straps, lambda, ...
        yresp,x,rpos,transformation,opt, ...
        y,rind,sind,normalise);
    
    lower_ci = confidencey(1,:,:);
    upper_ci = confidencey(2,:,:);
    
    % Check coverage
    for h = 1:H
        coverage(i,h) = true_irf(h) >= lower_ci(:,h) & true_irf(h) <= upper_ci(:,h);
    end
    temp_coverage = mean(coverage(1:i,:), 1);
    
    disp(['Current coverage is: ', num2str(temp_coverage)])
    % Calculate and display timing for this iteration
    timePerLoop(i) = toc(iterStart);
    estimatedTimeRemaining = mean(timePerLoop(1:i)) * (n_MC - i);
    
    fprintf('Iteration %d/%d completed in %.2f seconds\n', i, n_MC, timePerLoop(i));
    fprintf('Estimated time remaining: %.2f minutes\n', estimatedTimeRemaining/60);
end

% Results
mean_coverage = mean(coverage,1); % 1 = rows, to "collapse" the rows

fprintf('Empirical coverage rate: %.3f ', mean_coverage);

filename = "Bootstrap/results/0.5lambda_rind2_rzDGPVAR2_SVARLP_MC" + num2str(n_MC) + "_BS" + num2str(n_straps) + ".mat";

save(filename, "mean_coverage")

%% PLOT FOR PAPER - All Three
ci = 0.95;

plot_files = {
    "Bootstrap/results/newweightsT300s_rind2_rzDGPVAR2_SVARLP_MC1000_BS400.mat",
    "Bootstrap/results/rzDGPVAR2_SVAR_MC1000_BS400.mat",
    };

for i=1:length(plot_files)
    eval(['res' num2str(i) '= load(plot_files{i});'])
    eval(['coverage_data' num2str(i) '=res' num2str(i) '.mean_coverage;'])
end

fig = figure;
hold on
colors = ['g', 'b', 'k', 'm', 'c', 'y'];  % Basic MATLAB colors

plot(1:H, ci*ones(H,1), 'LineStyle','--', 'LineWidth', 1, 'Color','r');
for i=1:length(plot_files)
    eval(['coverage_data = coverage_data' num2str(i)]); 
    plot(1:H, coverage_data, 'LineStyle', '-', 'LineWidth', 1, 'Color', colors(i));
end
legend('95%', 'MC Coverage')
xlim([0,H])
ylim([0,1])
xlabel('Horizon of IRF')
ylabel('Average Coverage')
title('Monte Carlo Simulation of Model Average CIs')
legend('95%', 'Model Average Coverage', 'SVAR Coverage')

n_MC = 1000;
n_straps = 400;
subtitle(['Number of MC sims: ', num2str(n_MC), '; Number of bootstrap reps: ', num2str(n_straps)]);

saveas(fig, 'Fig/all_estimator_coverage_rzgdp_nov6.png')

%% DEBUGGING FROM NOW ON
% (Try a single MC run DGP simulation)

for i = 1:1
    close all
    T_total = 400;
    T_burn = 100;
    T_sample = T_total - T_burn;
    nlag = optlag;
    bootstrap_ci_method = "normal";
    clevel = 1.96; % we are interested 95% CI
    lambda = 0.5;
    timePerLoop = zeros(n_MC, 1);  % Store time for each iteration
    n_straps = 500;
    n_MC = 300;
    
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
    
    use_bootstrap = 1;
    normalise = 1;
    y = Y;
    
    [linysvar, confidenceysvar] = linSVAR(y,H,nlag,rind,sind,clevel,use_bootstrap,n_straps,normalise, bootstrap_ci_method);
    
    %[
    
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
    lp_use_boot = 0;
    
    [linylp, confidenceylp, ~, ~] = linlp(yresp, x, H, rpos, transformation, clevel, opt, lp_use_boot,nlag, n_straps, 'normal');
    
    disp(['svar h=0 response ', num2str(linysvar(1))])
    disp(['lp h=0 response ', num2str(linylp(1))])
    
    close all
    yfig = figure;
    hold on
    plot(Y(:,1))
    plot(Y(:,2))
    plot(Y(:,3))

    legend('Y_1', 'Y_2', 'Y_3')
    title('Sample data from Monte Carlo process using calibrated DGP')
%     plot(linysvar, 'Marker', 'X')
%     % plot(linylp, 'Marker','o')
%     plot(squeeze(confidenceysvar(1,:,:)))
%     plot(squeeze(confidenceysvar(2,:,:)))
%     legend('True IRF', 'SVAR')
    pause(2)
end

saveas(yfig, 'Fig/calibrated_dgp_ysample.png')


