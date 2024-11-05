
rng(123);

% Simulation to test the optimal forecast weights
% DGP for bivariate VAR(1): BY_t = A Y_{t-1} + epsilon
B = [1 0; 0.5 3];
A = [0.9 0; 0.5 0.5];
omega = eye(2); % variance of epsilon

% IRF parameters

% Simulation parameters
T = 200;              % Sample size
H = 20;               % forecast horizon
n_vars = 2;           % Number of variables
rind = 2;             % response variable index
sind = 1;             % shock variable index

Binv = inv(B);
sigma = Binv * omega * Binv';  % Reduced form covariance matrix

% Create the dataset
Y = zeros(T,n_vars);
eps = mvnrnd(zeros(n_vars,1), omega, T);

for t = 2:T
    By_t = Y(t,:) + Y(t-1,:) * A' + eps(t,:);
    Y(t,:) = By_t * inv(B);
end

residuals = zeros(T,n_vars);
for t = 2:T
    Yhat = Y(t-1,:) * A';
    residuals(t,:) = Y(t) - Yhat;
end
epsilon = residuals * inv(B)';

sep_idx = ceil(0.9*T);
train_y = Y(1:end-H,:);
test_y = Y(end-H+1:end,:);

% first and second inputs to this must have product H
W = generateWeightMatrix_v2(10,2,0.25);

% Forecasting and eval
cast_var = rind;
est_nlag = 1; % the number of lags used in estimation / forecasting

lp_forecasts = computeLPForecasts(train_y, est_nlag, H, cast_var);
svar_forecasts = computeSVARForecasts(train_y, est_nlag, H, cast_var);

msfe_matrix = zeros(1, size(W,2));
for i = 1:size(W,2)
    disp(['Currently testing weight set ', num2str(i), ' out of ', num2str(size(W,2))])
    weights = W(:,i);
    avg_forecast = weights.*svar_forecasts + (1-weights).*lp_forecasts;
    target_y = test_y(:,cast_var);
    msfe = (1/H) * sum((avg_forecast - target_y).^2);
    msfe_matrix(:,i) = msfe;
end

% Find the optimal weights
[min_msfe, opt_idx] = min(msfe_matrix);
optimal_weights = W(:, opt_idx);

%%
save('optimal_weights_block10size2_inc0.25', 'optimal_weights')

%%
bar(optimal_weights)
title('Weight on SVAR')
