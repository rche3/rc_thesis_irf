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
    x_la = [xreg(2:end,:), xorig_lagaug(2:end,:)];
    rpos=2;
end

linres = struct(); % struct to store linear model results

% compute the IRFs now with various methods
irfTypes = {'LP', 'LP_BC', 'LP_Penalised', 'LP_Lagaug', 'SVAR', 'VARLP_Avg'};
irfTypes = {'LP', 'SVAR', 'SVARLPAvg'};

L = length(irfTypes);

%% LP Estimation w/ DW Bootstrap
bootstrap=1;
nstraps = 100;
method = 'normal';
lambda = 0.5;
p=nlag;

%%% LP
e=1; est=irfTypes{e};
[linres.IRF.(est), linres.CI.(est)]=linlp(data,x,hor,rpos,transformation, clevel, opt, bootstrap, p, nstraps, method); 
%% VAR Estimation w/ Bootstrap
e=2; est=irfTypes{e};

% var preliminaries
normalise = 1;
y = [x(:,rpos), data];
p=nlag;
rind = [2,3]; % location of response variables in y
sind = 1; % location of shock variable in y
nstraps = 1000;
bootstrap = 1;
[linres.IRF.(est), linres.CI.(est), ~] = linSVAR(y, hor, p, rind, sind, clevel, bootstrap, nstraps, normalise);

%% VAR / LP Averaged
e=3; est = irfTypes{e};

% varlp setup
method = 'normal';
lambda = 0.5; % weight for VAR irfs
nstraps = 100;
[linres.IRF.(est), linres.CI.(est)] = linSVARLP_avg(hor,p, method, clevel, nstraps, lambda, ...
    data,x,rpos,transformation,opt, ...
    y, rind, sind, normalise);

%%
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
        plot(1:1:hor, zz, 'k-', 'HandleVisibility','off')
        hold on
        plot(1:1:hor, linres.IRF.(est)(j,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Point IRF');
        plot(1:1:hor, squeeze(linres.CI.(est)(1,j,:)), 'b--', 'DisplayName', '95% CI'); % lower bound CI
        plot(1:1:hor, squeeze(linres.CI.(est)(2,j,:)), 'b--', 'HandleVisibility', 'off'); % upper bound CI
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

sgtitle('Linear IRF Estimates')
set(linfig, 'Position', [100, 100, 1200, 800])
% saveas(linfig,'fig/linear_varlpavg_est.png')

