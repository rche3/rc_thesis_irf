%%% SCRIPT DESIGNED TO COMPUTE AVG/LP/SVAR Multipliers
clear all
close all
%clc
warning off;

addpath(genpath('Plots'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT SPECIFICATION CHOICES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timesample=1; % 1 means 1889-2015; 2 means post WWII
statechoice=2;  % 1 means default unemp threshold; 2 means default ZLB threshold
shockchoice=1; % 1 means news shock; 2 means BP shock
transformation=1; % 1 means Gordon-Krenn; 2 means Hall-Barro-Redlick
taxchoice=0; %0 means no taxes as control; 1 means ad taxy as a control
datafirstdiff = 0; % to difference the GDP / Gov Spending Data

trend=0; %4 means quartic; 2 means quadratic; 0 means no trend  
nlag=4;% number of lags of control variables
opt=0; % 0 means Newey-West with bandwidth=i; 1 means optimal bandwidth (takes much longer to run)'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if transformation==1
    disp('Gordon-Krenn')
elseif transformation==2
    disp('Hall-Barro-Redlick')
end

if shockchoice==1
    disp('News shock')
elseif shockchoice==2
    disp('BP shock')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hor=20; % horizon for IRFs
clevel= 1.96; % confidence level % 1.96  1.64;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run('figure_irfs_multipliers_setup.m')

% SETUP DATA
if trend==4
    xreg = [constant, t', tsq', tcu', tqu', shock];
    x = [xreg,  xorig]; %complete RHS
    x_la = [xreg(2:end, :), xorig_lagaug(2:end,:)];
    rpos=6; % position of shock
elseif trend==2
    xreg = [constant, t', tsq', shock];
    x= [xreg, xorig];
    x_la = [xreg(2:end, :), xorig_lagaug(2:end,:)];
    rpos=4;
elseif trend==0
    xreg = [constant, shock];
    x=[xreg,  xorig];
    x_la = [xreg(2:end,:), xorig_lagaug(2:end,:)];
    rpos=2;
end

%% Do some optimisation
obsshocksim = readtable('ObsShock_loss_bias_regcat.csv');

obsshocksim = obsshocksim(11:end,:);

% Define the error vectors
obsshocksim_bias_lp = obsshocksim.lp(1:end-1);
obsshocksim_bias_svar = obsshocksim.svar(1:end-1);

obsshocksim_bias_lp = abs(obsshocksim_bias_lp);
obsshocksim_bias_svar = abs(obsshocksim_bias_svar);

%%% OPTIMISE
options = optimoptions('fmincon', 'Display', 'iter');
objective_wrapper = @(x) objective(x, obsshocksim_bias_lp, obsshocksim_bias_svar);
[optimal_lambdas, fval] = fmincon(objective_wrapper, [0.5, 0.5, 0.5, 0.5], [], [], [], [], [0, 0, 0, 0], [1, 1, 1, 1], [], options);

% Calculate the weighted sum with optimal lambdas
optimal_weighted_sum = zeros(size(obsshocksim_bias_lp));
for i = 1:4
    idx = (1:5) + (i-1)*5;
    optimal_weighted_sum(idx) = optimal_lambdas(i) * obsshocksim_bias_lp(idx) + (1-optimal_lambdas(i)) * obsshocksim_bias_svar(idx);
end

% Print results
fprintf('Optimal lambdas:\n');
for i = 1:4
    fprintf('Lambda %d: %.4f\n', i, optimal_lambdas(i));
end
fprintf('Minimum weighted sum norm: %.4f\n', fval);

% Plot results
figure;
bar(1:4, optimal_lambdas);
xlabel('Segment');
ylabel('Optimal Lambda');
title('Optimal Lambdas for Each Segment');
ylim([0, 1]);

% Plot original vectors and optimal weighted sum
figure;
plot(obsshocksim_bias_lp, 'b-', 'LineWidth', 2);
hold on;
plot(obsshocksim_bias_svar, 'r-', 'LineWidth', 2);
plot(optimal_weighted_sum, 'g-', 'LineWidth', 2);
xlabel('Index');
ylabel('Value');
title('Original Vectors and Optimal Weighted Sum');
legend('obsshocksim\_bias\_lp', 'obsshocksim\_bias\_svar', 'Optimal Weighted Sum');
grid on;

% Add vertical lines to separate segments
for i = 1:3
    xline(i*5, '--k');
end

%% Linear LP, SVAR, Avg Point estimates

close all
linres = struct(); % struct to store linear model results
irfTypes = {'LP', 'SVAR', 'SVARLPAvg'};
L = length(irfTypes);

%%% do estimation
% (LP)
e=1; est=irfTypes{e};
bootstrap=0;
nstraps = 100;
method = 'normal';
lambda = 0.5;
p=nlag;
[linres.IRF.(est), linres.CI.(est)]=linlp(data,x,hor,rpos,transformation, clevel, opt, bootstrap, p, nstraps, method); 
% (SVAR)
e=2; est=irfTypes{e};
normalise = 1;
y = [x(:,rpos), data];
p=nlag;
rind = [2,3]; % location of response variables in y
sind = 1; % location of shock variable in y
nstraps = 1000;
bootstrap = 1;
[linres.IRF.(est), linres.CI.(est), ~] = linSVAR(y, hor, p, rind, sind, clevel, bootstrap, nstraps, normalise);
% (Model average)
e=3; est=irfTypes{e};
horizon_weights = [0.5 * ones(1, hor)]; % weighting on LP
linres.IRF.(est) = horizon_weights.*linres.IRF.LP + (1-horizon_weights).*linres.IRF.SVAR;

% convert to multiplier
for i=1:L
    est = irfTypes{i};
    linres.cum_mult.(est) = cumsum(linres.IRF.(est)(2,:))./cumsum(linres.IRF.(est)(1,:));
end

run('setup.m')

ax = zeros(1, length(irfTypes));

fig = figure;
for i = 1:length(irfTypes)
    irfType = irfTypes{i};
    ax(i) = subplot(1,3,i);
    hold on 
    plot(1:1:hor, linres.cum_mult.(irfType), 'Color', line_colors{i}, 'LineWidth', 1.5, 'DisplayName', irfType); % plot point estimates for multiplier
    axis tight
    legend('Location', 'east')
end

linkaxes(ax, 'xy');  % This links both x and y axes

sgtitle('Linear model: cumulative spending multiplier (point estimates) with different estimators');

fig.Position = [100 100 800 300];  % [left bottom width height]

saveas(fig, 'multiplier_point_est.png')

disp(['Long run LP multiplier is: ' num2str(linres.cum_mult.LP(end))])
disp(['Long run SVAR multiplier is: ' num2str(linres.cum_mult.SVAR(end))])
disp(['Long run Model Averaged multiplier is: ' num2str(linres.cum_mult.SVARLPAvg(end))])

%%
lpmult = linres.cum_mult.LP(end);
svarmult = linres.cum_mult.SVAR(end);
mdlavgmult = linres.cum_mult.SVARLPAvg(end);

[lpmult; svarmult; mdlavgmult]

%% State dependent LP point estimates
z=xorig; % lagged observations
In=fu(nlag+1:end); % truncated indicator variable for LP regression
shock = shock; % from the setup script
bootstrap = 0;
nstraps = 200;
method = 'normal';

% First, normal regression SEs
[rzirfa, rzirfb, rzcia, rzcib, ~, ~] = statelp(data,shock,z,In,hor, ...
    transformation,clevel,opt, ...
    bootstrap,method,nlag,nstraps);

cum_mult_statea= cumsum(stateay(2,:))./cumsum(stateay(1,:));
cum_mult_stateb= cumsum(stateby(2,:))./cumsum(stateby(1,:));

%%
function weighted_sum_norm = objective(lambdas, vec1, vec2)
    weighted_sum = zeros(size(vec1));
    for i = 1:4
        idx = (1:5) + (i-1)*5;
        weighted_sum(idx) = lambdas(i) * vec1(idx) + (1-lambdas(i)) * vec2(idx);
    end
    weighted_sum_norm = norm(weighted_sum);
end