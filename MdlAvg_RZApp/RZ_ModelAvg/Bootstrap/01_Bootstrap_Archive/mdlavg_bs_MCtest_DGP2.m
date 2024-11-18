clc; clear; close all;

%%%%%
% Script to use Monte Carlo simulation to test the validity of the dependent 
% wild bootstrapping for the LP and SVAR model average IRF estimator

% Use a VAR(4) model as the DGP
%%%%%

% Set random seed for reproducibility
rng(123);

% Simulation parameters
T = 200;              % Sample size
H = 20;               % IRF horizon
n_vars = 3;           % Number of variables
rind = 3;             % response variable index
sind = 1;             % shock variable index

% True STRUCTURAL parameters
B = [1, 0, 0;         % variable 1 can only be contemporaneously shocked by itself
     -0.3, 1, 0;      % variable 2 can be shocked by itself + variable 1
     0.2, -0.1, 1];   % variable 3 can be shocked by itself + all other variables

% True structural VAR coefficients
A1 = [0.30, 0.15, -0.08;    
      0.1, 0.26, 0.15;      
     -0.08, 0.1, 0.25];     

A2 = [0.20, 0.08, -0.04;    
      0.04, 0.20, 0.08;     
     -0.04, 0.04, 0.20];    

A3 = [0.10, 0.04, -0.02;   
      0.02, 0.10, 0.04;    
     -0.02, 0.02, 0.10];  

A4 = [0.05, 0.02, -0.01;  
      0.01, 0.05, 0.02;   
     -0.01, 0.01, 0.05];  

% Structural shock variances (diagonal matrix)
D = eye(3);

% Construct companion matrix for VAR(4)
A = [A1, A2, A3, A4;
    eye(3), zeros(3,9);
    zeros(3), eye(3), zeros(3,6);
    zeros(3), zeros(3), eye(3), zeros(3,3)];

disp(['Max eigenvalue is ' num2str(max(abs(eig(A))))])


% Structural shock variances (diagonal matrix)
D = eye(3);  % Different variances for structural shocks - just set them to be unit variance

% Implied reduced form parameters
Binv = inv(B);
Sigma = Binv * D * Binv';  % Reduced form covariance matrix

% Compute true structural IRFs
true_irf = zeros(H, 1);
impact = Binv(:,sind); % Initial impact of structural shock sind
true_irf(1) = impact(rind);

% Compute subsequent horizons using companion matrix
state = [impact; zeros(9,1)]; % State vector includes 4 lags of 3 variables
for h = 2:H
    state = A * state;
    true_irf(h) = state(rind);
end

%% Monte Carlo simulation - LP
% specs for model average estimation
nlag = 4;
bootstrap_ci_method = "normal";
use_bootstrap = 1;
clevel = 1.96; % we are interested 95% CI
n_straps = 200;
n_MC = 500;
timePerLoop = zeros(n_MC, 1);  % Store time for each iteration
coverage = zeros(n_MC, H);

for i = 1:n_MC
    % Generate data from VAR(1)
    iterStart = tic;  % Start timing this iteration

    disp(['Running MC rep number ', num2str(i)]);
    Y = zeros(T, 3);
    eps = mvnrnd(zeros(3,1), D, T);

    for t = 5:T
        By_t = eps(t,:);
        for lag = 1:4
            By_t = By_t + Y(t-lag,:) * eval(['A' num2str(lag)])';
        end
        Y(t,:) = By_t * inv(B);
    end
    

    % get structural shocks for LP estimation
    epsilon = nan(T, n_vars);
    for t = 5:T
        pred = zeros(1, n_vars);
        for lag = 1:4
            pred = pred + Y(t-lag,:) * eval(['A' num2str(lag)])';
        end
        resids = Y(t,:) - pred;
        epsilon(t,:) = resids * B';
    end

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

filename = "Bootstrap/results/VAR4_lp_MC" + num2str(n_MC) + "_BS" + num2str(n_straps) + ".mat";
save(filename, "mean_coverage")

%% Monte Carlo simulation - SVAR with bootstrap

% specs for model average estimation
nlag = 4;
bootstrap_ci_method = "normal";
clevel = 1.96; % we are interested 95% CI
normalise = 1;
use_bootstrap = 1;
n_straps = 200;
n_MC = 500;
timePerLoop = zeros(n_MC, 1);  % Store time for each iteration
coverage = zeros(n_MC, H);

for i = 1:n_MC
    % Generate data from VAR(1)
    iterStart = tic;  % Start timing this iteration

    disp(['Running MC rep number ', num2str(i)]);
    Y = zeros(T, 3);
    eps = mvnrnd(zeros(3,1), D, T);
    for t = 5:T
        By_t = eps(t,:);
        for lag = 1:4
            By_t = By_t + Y(t-lag,:) * eval(['A' num2str(lag)])';
        end
        Y(t,:) = By_t * inv(B);
    end

    % VAR specs
    y = Y;
    
    % Estimate linSVARLP_avg on MC sample
    [~, confidencey] = linSVAR(y,H,nlag,rind,sind,clevel,use_bootstrap,n_straps,normalise);
    
    lower_ci = confidencey(1,:,:);
    upper_ci = confidencey(2,:,:);
    
    % Check coverage
    for h = 1:H
%         if h < 5
%             disp([liny(h), true_irf(h)])
%         else
%         end
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

filename = "Bootstrap/results/VAR4_SVAR_MC" + num2str(n_MC) + "_BS" + num2str(n_straps) + ".mat";

save(filename, "mean_coverage")

%% Monte Carlo simulation - Model average with bootstrap
% specs for model average estimation
nlag = 4;
bootstrap_ci_method = "normal";
use_bootstrap = 1;
clevel = 1.96; % we are interested 95% CI
n_straps = 300;
n_MC = 500;
lambda = 0.5;
timePerLoop = zeros(n_MC, 1);  % Store time for each iteration
coverage = zeros(n_MC, H);

for i = 1:n_MC
    % Generate data from VAR(1)
    iterStart = tic;  % Start timing this iteration

    disp(['Running MC rep number ', num2str(i)]);
    Y = zeros(T, 3);
    eps = mvnrnd(zeros(3,1), D, T);

    for t = 5:T
        By_t = eps(t,:);
        for lag = 1:4
            By_t = By_t + Y(t-lag,:) * eval(['A' num2str(lag)])';
        end
        Y(t,:) = By_t * inv(B);
    end

    % get structural shocks for LP estimation
    epsilon = nan(T, n_vars);
    for t = 5:T
        pred = zeros(1, n_vars);
        for lag = 1:4
            pred = pred + Y(t-lag,:) * eval(['A' num2str(lag)])';
        end
        resids = Y(t,:) - pred;
        epsilon(t,:) = resids * B';
    end

    % LP
    yresp = Y(:,rind);
    shock = epsilon(:,sind);
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
    y2y3 = Y(:, [rind-1, rind]);
    y2y3 = y2y3(nlag+1:end,:);
    y1 = shock(nlag+1:end,:);
    y = [y1, y2y3];
    
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

filename = "Bootstrap/results/VAR4_mdlavg_MC" + num2str(n_MC) + "_BS" + num2str(n_straps) + ".mat";
save(filename, "mean_coverage")

%% Plotting coverage rate
ci = 0.95;

res = load("Bootstrap/results/VAR4_mdlavg_MC500_BS300.mat");
coverage_data = res.mean_coverage

fig = figure;
hold on
plot(1:H, ci*ones(H,1), 'LineStyle','--', 'LineWidth', 1, 'Color','r');
plot(1:H, coverage_data, 'LineStyle', '-', 'LineWidth', 1, 'Color','black')
legend('95%', 'MC Coverage')
xlim([0,H])
ylim([0,1])
xlabel('Horizon of IRF')
ylabel('Average Coverage')
title('Monte Carlo Simulation of Model Averaged CIs')

n_MC = 500;
n_straps = 300;

subtitle(['Number of MC sims: ', num2str(n_MC), '; Number of bootstrap reps: ', num2str(n_straps)]);

saveas(fig, 'Fig/MC_test_mdlavgbootstrap_var4.png')

%% Debug
Y = zeros(T, 3);
eps = mvnrnd(zeros(3,1), D, T);
for t = 5:T
    By_t = eps(t,:);
    for lag = 1:4
        By_t = By_t + Y(t-lag,:) * eval(['A' num2str(lag)])';
    end
    Y(t,:) = By_t * inv(B);
end
y = Y;

epsilon = nan(T, n_vars);
for t = 5:T
    pred = zeros(1, n_vars);
    for lag = 1:4
        pred = pred + Y(t-lag,:) * eval(['A' num2str(lag)])';
    end
    resids = Y(t,:) - pred;
    epsilon(t,:) = resids * B';
end

% SVAR
nlag = 4;
bootstrap_ci_method = "normal";
clevel = 1.96; % we are interested 95% CI
normalise = 1;
n_straps = 200;
use_bootstrap = 1;
n_MC = 200;
timePerLoop = zeros(n_MC, 1);  % Store time for each iteration
coverage = zeros(n_MC, H);

[linysvar, confidenceysvar] = linSVAR(y,H,nlag,rind,sind,clevel,use_bootstrap,n_straps,normalise);

var4model = varm(3, 4);
estvar4model = estimate(var4model,y);
res = irf(estvar4model);
linysvar = res(:, sind, rind);

 % LP
yresp = Y(:,rind);
shock = epsilon(:,sind);
ylag = lagmatrix(Y, [1:nlag]);
ylag = ylag(nlag+1:end,:);
Teff = size(ylag, 1);

x = [ones(Teff,1) shock(nlag+1:end,:) ylag];
yresp = yresp(nlag+1:end,:);

rpos = 2;
transformation = 1;
opt = 0;

[linylp, confidenceylp] = linlp(yresp, x, H, rpos, transformation, clevel, opt, 0, nlag, n_straps, bootstrap_ci_method);

figure;
hold on
plot(true_irf)
plot(linysvar, 'Marker', 'X')
plot(confidenceysvar(1,:))
plot(confidenceysvar(2,:))
% plot(linylp, 'Marker', 'X')
% plot(confidenceylp(1,:))
% plot(confidenceylp(2,:))
legend('True IRF', 'IRF Point Est SVAR','SVAR Lower', 'SVAR Upper', 'IRF Point Est LP', 'LP Lower', 'LP Upper')


