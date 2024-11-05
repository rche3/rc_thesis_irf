clc; clear; close all;

%%%%%
% Script to use Monte Carlo simulation to test the validity of the dependent 
% wild bootstrapping for the LP and SVAR model average IRF estimator

% Use a VAR(1) model as the DGP
%%%%%

% Set random seed for reproducibility
rng(123);

% Simulation parameters
T_total = 400;              % Sample size
T_burn = 100;
H = 20;               % IRF horizon
n_vars = 3;           % Number of variables
rind = 3;             % response variable index
sind = 1;             % shock variable index

% True STRUCTURAL parameters
B = [1, 0, 0;         % variable 1 can only be contemporaneously shocked by itself
     -0.3, 1, 0;      % variable 2 can be shocked by itself + variable 1
     0.2, -0.1, 1];   % variable 3 can be shocked by itself + all other variables

% True structural VAR coefficients
A_str = [0.5, 0.2, -0.1;
         0.1, 0.4, 0.2;
         -0.1, 0.1, 0.6];

% Structural shock variances (diagonal matrix)
D = eye(3);  % Different variances for structural shocks - just set them to be unit variance

% Implied reduced form parameters
A = A_str;  % In this case, coefficients are same (you could make them different)
Binv = inv(B);
Sigma = Binv * D * Binv';  % Reduced form covariance matrix

% Compute true structural IRFs
true_irf = zeros(H, 1);
impact = Binv(:,sind);  % Initial impact of structural shock sind
true_irf(1) = impact(rind);  % Store response of variable rind
companion = A;

% Compute subsequent horizons for true IRF
for i = 2:H
    temp = companion^(i-1) * impact;
    true_irf(i) = temp(rind);
end

%% Monte Carlo simulation - LP
% specs for model average estimation
nlag = 1;
bootstrap_ci_method = "normal";
use_bootstrap = 0;
clevel = 1.96; % we are interested 95% CI
n_straps = 200;
n_MC = 3000;
timePerLoop = zeros(n_MC, 1);  % Store time for each iteration

for i = 1:n_MC
    % Generate data from VAR(1)
    iterStart = tic;  % Start timing this iteration

    disp(['Running MC rep number ', num2str(i)]);
    Y_total = zeros(T_total, 3);
    eps = mvnrnd(zeros(3,1), D, T_total);
    for t = 2:T_total
        By_t = Y_total(t,:) + Y_total(t-1,:) * A' + eps(t,:);
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

filename = "Bootstrap/results/newversionVAR1_lp_MC" + num2str(n_MC) + "_BS" + num2str(n_straps) + ".mat";
save(filename, "mean_coverage")

%% Monte Carlo simulation - SVAR with bootstrap

% specs for model average estimation
nlag = 1;
bootstrap_ci_method = "normal";
clevel = 1.96; % we are interested 95% CI
timePerLoop = zeros(n_MC, 1);  % Store time for each iteration
n_straps = 200;
n_MC = 300;

for i = 1:n_MC
    % Generate data from VAR(1)
    iterStart = tic;  % Start timing this iteration

    disp(['Running MC rep number ', num2str(i)]);
    Y_total = zeros(T_total, 3);
    eps = mvnrnd(zeros(3,1), D, T_total);
    for t = 2:T_total
        By_t = Y_total(t,:) + Y_total(t-1,:) * A' + eps(t,:);
        Y_total(t,:) = By_t * inv(B);
    end
    Y = Y_total(T_burn+1:end,:);

    % prepare variables for estimation
    data = Y(:,rind);
    % not sure how to convert this here given we are no longer reg on a shock
    
    shock = Y(:,sind); % the shock here should be initialised as contemporaneous to the response variable
    ylag = lagmatrix(Y, [1:nlag]);
    ylag = ylag(nlag+1:end,:);
    Teff = size(ylag, 1);

    % create final Y and X data vectors for LP
    rpos = rind;    
    use_bootstrap = 1;
    
    % Estimate linlp on MC sample
    % VAR specs
    normalise = 1;
    y = Y;
    
    % Estimate linSVARLP_avg on MC sample
    [liny, confidencey] = linSVAR(y,H,nlag,rind,sind,clevel,use_bootstrap,n_straps,normalise, bootstrap_ci_method);
    
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
%     pause(5)
end

% Results
mean_coverage = mean(coverage,1); % 1 = rows, to "collapse" the rows

filename = "Bootstrap/results/newversionVAR1_SVAR_MC" + num2str(n_MC) + "_BS" + num2str(n_straps) + ".mat";

save(filename, "mean_coverage")

%% Monte Carlo simulation - model average
% settings
nlag = 1;
bootstrap_ci_method = "normal";
clevel = 1.96; % we are interested 95% CI
lambda = 0.5;
timePerLoop = zeros(n_MC, 1);  % Store time for each iteration
n_MC = 500;
n_straps = 200;

for i = 1:n_MC
    iterStart = tic;
    disp(['Running MC rep number ', num2str(i)]);
    % Generate data from VAR(1)
    Y_total = zeros(T_total, 3);
    eps = mvnrnd(zeros(3,1), D, T_total);
    for t = 2:T_total
        By_t = Y_total(t,:) + Y_total(t-1,:) * A' + eps(t,:);
        Y_total(t,:) = By_t * inv(B);
    end
    Y = Y_total(T_burn+1:end,:);
    
    % LP
    yresp = Y(:,rind);
    shock = Y(:,sind);
    ylag = lagmatrix(Y, [1:nlag]);
    ylag = ylag(nlag+1:end,:);
    Teff = size(ylag, 1);
    
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
%%

filename = "Bootstrap/results/newversionVAR1_SVARLP_MC" + num2str(n_MC) + "_BS" + num2str(n_straps) + ".mat";

save(filename, "mean_coverage")

%% PLOT FOR PAPER
ci = 0.95;

res = load("Bootstrap/results/newversionVAR1_SVARLP_MC500_BS200.mat");
coverage_data = res.mean_coverage;

fig = figure;
hold on
plot(1:H, ci*ones(H,1), 'LineStyle','--', 'LineWidth', 1, 'Color','r');
plot(1:H, coverage_data, 'LineStyle', '-', 'LineWidth', 1, 'Color','black')
legend('95%', 'MC Coverage')
xlim([0,H])
ylim([0,1])
xlabel('Horizon of IRF')
ylabel('Average Coverage')
title('Monte Carlo Simulation of Model Average CIs')

n_MC = 500;
n_straps = 200;
subtitle(['Number of MC sims: ', num2str(n_MC), '; Number of bootstrap reps: ', num2str(n_straps)]);

saveas(fig, 'Fig/newversion_MC_test_mdlavg_strap.png')

%% Plotting point estimate to see if roughly correct
% (Try a single MC run DGP simulation)
nlag = 1;
bootstrap_ci_method = "normal";
clevel = 1.96; % we are interested 95% CI
lambda = 0.5;
timePerLoop = zeros(n_MC, 1);  % Store time for each iteration
n_straps = 200;
n_MC = 300;

Y = zeros(T, 3);
eps = mvnrnd(zeros(3,1), D, T);
for t = 2:T
    By_t = Y(t,:) + Y(t-1,:) * A' + eps(t,:);
    Y(t,:) = By_t * inv(B);
end

% get structural shocks for LP estimation
resids = zeros(T,3);
for t = 2:T
    Yhat = Y(t-1,:) * A';
    resids(t,:) = Y(t) - Yhat;
end
epsilon = resids * inv(B)';

% VAR 
use_bootstrap = 0;
normalise = 1;
y = Y;

[linysvar, confidenceysvar] = linSVAR(y,H,nlag,rind,sind,clevel,use_bootstrap,n_straps,normalise);

% LP
yresp = Y(:,rind);
shock = epsilon(:,sind);
ylag = lagmatrix(Y, [1:nlag]);
ylag = ylag(nlag+1:end,:);
Teff = size(ylag, 1);

% create final Y and X data vectors for LP
x = [ones(Teff,1) shock(nlag+1:end,:) ylag];
yresp = yresp(nlag+1:end,:);
rpos = 2;
transformation = 1;
opt = 0;

[linylp, confidenceylp, ~, ~] = linlp(yresp, x, H, rpos, transformation, clevel, opt, use_bootstrap,nlag, n_straps, 'normal');

disp(['svar h=0 response', num2str(linysvar(1))])
disp(['lp h=0 response', num2str(linylp(1))])

close all
figure;
hold on
plot(true_irf)
plot(linysvar, 'Marker', 'X')
plot(linylp, 'Marker','o')
legend('True IRF', 'Last IRF Point Est SVAR', 'Last IRF Point Est LP')