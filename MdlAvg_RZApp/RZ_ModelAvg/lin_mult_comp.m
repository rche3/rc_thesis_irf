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
hor=21; % horizon for IRFs
clevel= 1.96; % confidence level % 1.96  1.64;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run('figure_irfs_multipliers_setup.m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LINEAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    x_la = [xreg(2:end,:), xorig_lagaug(2:end,:)]; % this is for lag augmentation
    rpos=2;
end

linres = struct(); % struct to store linear model results

% compute the IRFs now with various methods
irfTypes = {'LP', 'SVAR', 'SVARLPAvg_50', 'SVARLPAvg_Opt'};
L = length(irfTypes);

bootstrap=1;
nstraps = 400;
method = 'percentile';
p=nlag;

%%% LP
e=1; est=irfTypes{e};
[linres.IRF.(est), linres.CI.(est), ~, ~, linres.mult_lpu.(est)]=linlp(data,x,hor,rpos,transformation, clevel, opt, bootstrap, p, nstraps, method); 

%%% VAR Estimation w/ Bootstrap
e=2; est=irfTypes{e};

normalise = 1;
y = [x(:,rpos), data];
rind = [2,3]; % location of response variables in y
sind = 1; % location of shock variable in y
[linres.IRF.(est), linres.CI.(est), ~, linres.mult_lpu.(est)] = linSVAR(y, hor, p, rind, sind, clevel, bootstrap, nstraps, normalise, method);

%%% VAR / LP Averaged w/ Custom Bootstrap
e=3; est = irfTypes{e};

lambda = 0.5*ones(1,hor); % weight for VAR irfs
[linres.IRF.(est), linres.CI.(est), linres.mult_lpu.(est)] = linSVARLP_avg(hor,p, method, clevel, nstraps, lambda, bootstrap, ...
    data,x,rpos,transformation,opt, ...
    y, rind, sind, normalise);

%%% VAR / LP Averaged w/ Custom Bootstrap, Opt
e=4; est = irfTypes{e};

lambda_read = load("Weight Selection/opt_weight.mat");
lambda = [0.5, lambda_read.opt_weights];
[linres.IRF.(est), linres.CI.(est), linres.mult_lpu.(est)] = linSVARLP_avg(hor,p, method, clevel, nstraps, lambda, bootstrap, ...
    data,x,rpos,transformation,opt, ...
    y, rind, sind, normalise);

irfTypes_concat = '';

for l = 1:L
    irfTypes_concat = [irfTypes_concat irfTypes{l} '_'];
end

filename = "Results/lin_" + method + "_" + irfTypes_concat + "BS" + num2str(nstraps);
save(filename, "linres")

%% Load data to plot

res = load("Results/lin_percentile_LP_SVAR_SVARLPAvg_50_SVARLPAvg_Opt_BS400.mat");
linres = res.linres;

% long run estimates for multiplier
disp('Long run multiplier estimates')
point = [];
for l=1:L
    est = irfTypes{l};
    point = [point; linres.mult_lpu.(est)(3,end)];
end
disp(point) 

disp('Long run multiplier CI')
ci = [];
for l=1:L
    est = irfTypes{l};
    ci = [ci; linres.mult_lpu.(est)(1,end) linres.mult_lpu.(est)(2,end)];
end
disp(ci) 

%% IRF Plots

close all
zz=zeros(1,hor);
linfig = figure;
ylabels = {'Government Spending', 'GDP'};

% Create an array to store axes handles
ax = [];

for l=1:L
    est = irfTypes{l};
    for j=1:2
        subplot(2,L,l+(j-1)*L)
        plot(0:1:hor-1, zz, 'k-', 'HandleVisibility','off')
        hold on
        plot(0:1:hor-1, linres.IRF.(est)(j,1:hor), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Point IRF');
        plot(0:1:hor-1, squeeze(linres.CI.(est)(1,j,1:hor)), 'b--', 'DisplayName', '95% CI'); % lower bound CI
        plot(0:1:hor-1, squeeze(linres.CI.(est)(2,j,1:hor)), 'b--', 'HandleVisibility', 'off'); % upper bound CI
        axis tight
        ylabel(ylabels{j})
        xlabel('Horizon (Quarters)')
        title(est)
        legend('show')
        
        % Store the current axes handle
        ax = [ax gca];
    end
end

% Link y-axes of all subplots
linkaxes(ax, 'y');
sgtitle(['Linear IRF Estimates', ' | Number of bootstrap replications: ', num2str(nstraps), ' | CI Method: ', method])
set(linfig, 'Position', [300, 300, 1200, 800])
saveas(linfig,'fig/linear_varlpavg_est.png')

%% Multiplier Plots
close all

h_math_start = 1; % the mathematical horizon h at which we want to start plotting multipliers

zz=zeros(1,hor-h_math_start);
linfig = figure;
ylabels = {'Government Spending', 'GDP'};

% Create an array to store axes handles
ax = [];
titles = {'Local Projection', 'SVAR', 'SVAR LP 50/50', 'SVAR LP Opt Weight'};

for l=1:L
    est = irfTypes{l};
    subplot(1,L,l)
    plot(h_math_start:1:hor-1, zz, 'k-', 'HandleVisibility','off')
    hold on
    plot(h_math_start:1:hor-1, linres.mult_lpu.(est)(3,1+h_math_start:hor), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Multiple Point Estimate');
    plot(h_math_start:1:hor-1, linres.mult_lpu.(est)(1,1+h_math_start:hor), 'b--', 'DisplayName', '95% CI'); % lower bound CI
    plot(h_math_start:1:hor-1, linres.mult_lpu.(est)(2,1+h_math_start:hor), 'b--', 'HandleVisibility', 'off'); % upper bound CI
    axis tight
    ylabel('Government Spending Multiplier')
    xlabel('Horizon (Quarters)')
    title(titles{l})
    legend('show')
    
    % Store the current axes handle
    ax = [ax gca];
    xlim([h_math_start, hor-1])
    xticks(h_math_start:2:hor-1)
end

% Link y-axes of all subplots
linkaxes(ax, 'y');

sgtitle(['Linear multiplier Estimates', ' | Number of bootstrap replications: ', num2str(nstraps), ' | CI Method: ', method])
set(linfig, 'Position', [300, 300, 1000, 400])
saveas(linfig,'Fig/linest_multiplier_perc.png')
