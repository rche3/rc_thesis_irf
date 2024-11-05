
rng(0);

% Simulation parameters
T_total = 400;
T_burn = 100;

H = 20;               % IRF horizon
n_vars = 3;           % Number of variables
rind = 2;             % response variable index
sind = 1;             % shock variable index

% VAR(p) DGP
[B, AA, ~, optlag] = generate_DGP_rz(10,0);
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
disp(['Max eigenvalue is ' num2str(max(abs(eig(A))))])

Binv = inv(B);

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
T = size(Y,1);

sep_idx = ceil(0.9*T); % this is also the number of observations for training
T_train = sep_idx;
num_oos = T - sep_idx; % this is the number of pseudo-OOS observations

opt_weights = zeros(H,1);
est_nlag = 2; % the number of lags used in estimation / forecasting
cast_var = rind;

weight_set = 0:0.1:1;

horizon_weight_msfe = zeros(size(weight_set, 2), H);

for h=1:H
    temp_msfe_h = zeros(length(weight_set),1);
    for i = 1:length(weight_set)
        lambda = weight_set(i);
        disp(['Testing horizon ', num2str(h), ' and weight ', num2str(lambda)])
        % Forecasting and eval
        msfe_counter = 0;
        total_forecasts = T - T_train;
        for t = T_train+h+1:T-h
            train_final_idx = t-h;
            train_start_idx = train_final_idx - T_train +1;
            train_y = Y(train_start_idx:train_final_idx,:);
            test_y = Y(t+h,:);
            lp_forecast = forecastLP(train_y, est_nlag, h, cast_var);
            svar_forecast = forecastVAR(train_y, est_nlag, h, cast_var);
            avg_forecast = lambda.*svar_forecast + (1-lambda).*lp_forecast;
            target_y = test_y(:,cast_var);
            msfe = (avg_forecast - target_y).^2;
            msfe_counter = msfe_counter + msfe;
        end
        temp_msfe_h(i) = msfe_counter / total_forecasts;
    end
    horizon_weight_msfe(:,h) = temp_msfe_h;
end

opt_weight_idx = zeros(1,H);
for h = 1:H
    horizon_msfe_vec = horizon_weight_msfe(:,h);
    disp(horizon_msfe_vec)
    [~, idx] = min(horizon_msfe_vec);
    opt_weight_idx(h) = idx;
end

opt_weights = weight_set(opt_weight_idx);
opt_weights(1)=0.5;

disp(opt_weights')
bar(opt_weights)
title('Weight on SVAR')

save('Weight Selection/opt_weight.mat', "opt_weights")

