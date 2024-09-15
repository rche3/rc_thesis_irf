clear all
close all
%clc
warning off;

addpath(genpath('Plots'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT SPECIFICATION CHOICES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timesample=1; % 1 means 1889-2015; 2 means post WWII
statechoice=1;  % 1 means default unemp threshold; 2 means default ZLB threshold
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

%% test IRFs
y = [x(:,rpos), data];
p = 4;
[T, k] = size(y);

model = varm(k, p);
[mdl, ~, ~, res] = estimate(model, y);
[irf_var, lower, upper] = irf(mdl, E=res);

% lpw var script
[irf_var_lpw, ~] = linSVAR(data, x, hor, rpos, transformation, clevel, opt, nlag, 0);

% lp comparison
[irf_lp, ~] = linlp(data,x,hor,rpos,transformation, clevel, opt, nlag, 0);

% plot
figure;

var_irfshock1resp2 = irf_var(:,1,2)/irf_var(1,1,1);
lpwvar_irfshock1resp2 = irf_var_lpw(1,:);
lp_irfshock1resp2 = irf_lp(1,:);


figure;
hold on
subplot(1,3,1)
plot(0:19, var_irfshock1resp2)
subplot(1,3,2)
plot(0:19, lpwvar_irfshock1resp2)
subplot(1,3,3)
plot(0:19, lp_irfshock1resp2)
sgtitle('VAR MATLAB, Manual VAR - LPW, LP')


% for i = 1:3
%     subplot(3,1,i)
%     plot(1:20, squeeze(irf_data(i,1,:)), 'b-', ...
%          1:20, squeeze(lower(i,1,:)), 'r--', ...
%          1:20, squeeze(upper(i,1,:)), 'r--')
%     title(['Response of ', variable_names{i}, ' to Shock in Variable 1'])
%     xlabel('Periods')
%     ylabel('Response')
%     grid on
% end

%%
[irf_var_lpw, ~] = linSVAR(data, x, hor, rpos, transformation, clevel, opt, nlag, 0);
lpwvar_irfshock1resp2 = irf_var_lpw(1,:);

figure;
plot(0:19, lpwvar_irfshock1resp2)

%% 
clear all;
close all;
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
hor=20; % horizon for IRFs
clevel= 1.96; % confidence level % 1.96  1.64;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run('figure_irfs_multipliers_setup.m')

% compute the non-linear averaged estimator + for the model averaged forecast
if trend==4 % allow quartic trend
    x=[t', tsq', tcu', tqu', constant, (1-fu).*constant, fu.*shock, (1-fu).*shock, repmat(fu,1,size(x,2)).*x, repmat((1-fu),1,size(x,2)).*x];
    rpost=7; 
elseif trend==2 % allow quadratic trend
    x=[t', tsq', constant, (1-fu).*constant, fu.*shock, (1-fu).*shock, repmat(fu,1,size(x,2)).*x, repmat((1-fu),1,size(x,2)).*x];
    rpost=5; 
elseif trend==0 %no trend
    x=[constant, (1-fu).*constant, fu.*shock, (1-fu).*shock, repmat(fu,1,size(x,2)).*x, repmat((1-fu),1,size(x,2)).*x];
    rpost=3; % position of shock
end

%% separate TVAR and State Dependent LP
% FIRST, TEST OUR TVAR (SHOWED THIS CODE WORKS)
y = y_tvar; % ordered as shock, gov, gdp
I = fu;
p=4;
hor=20;
rind=[2,3]; % irfs to gov and GDP
sind=1; % shock ordered first in the y vector

[irfa, irfb, ~, ~, ~, ~, beta] = stateSVAR(y, I, p, hor,rind,sind);


%% model averaged estimator - point






[stateay, stateby, confidenceya, confidenceyb, bs_betadista, bs_betadistb, pvala, pvalb] = stateSVARLP_avg(data, ...
    x,hor,rpost,transformation, clevel, opt, bootstrap, emp, nlag, nstraps)


%% model averaged estimator - bootstrap / p val


%% test state dependency
close all
shockchoice=1; % 1 means news shock; 2 means BP shock

if shockchoice==1
    cifile_rec= xlsread('TVAR_output.xlsx',1);
    cifile_zlb= xlsread('TVAR_output.xlsx',2);
else
    cifile_rec= xlsread('TVAR_output.xlsx',3);
    cifile_zlb= xlsread('TVAR_output.xlsx',4);
end

if shockchoice==1
    recirf=cifile_rec(28:end, 2+3:10);
    nonrecirf=cifile_rec(28:end,13+3:21);
    
    zlbirf=cifile_zlb(28:end-1, 2+3:10);
    nonzlbirf=cifile_zlb(28:end-1,13+3:21);
else
    recirf=cifile_rec(28:end-1, 2:7);
    nonrecirf=cifile_rec(28:end-1,11:16);
    
    zlbirf=cifile_zlb(28:end-1, 2:7);
    nonzlbirf=cifile_zlb(28:end-1,10:15);
end

hor=20;

zz=zeros(20,1);

subplot(2,2,1)
hold on
plot(1:1:hor, zz, 'k-')
plot(1:1:hor, nonzlbirf(:,1), 'b--','LineWidth', 1.5)
plot(1:1:hor, zlbirf(:,1), 'r-o', 'LineWidth', 1.5);
xlabel('quarter')
title('ZLB-dependent: Government spending - RZ')
axis tight
legend('y=0', 'Non ZLB', 'ZLB')

subplot(2,2,2)
hold on
plot(1:1:hor, zz, 'k-')
plot(1:1:hor, irfa(1,:), 'b--','LineWidth', 1.5)
plot(1:1:hor, irfb(1,:), 'r-o', 'LineWidth', 1.5);
xlabel('quarter')
title('ZLB-dependent: Government spending - my TVAR ')
axis tight
legend('y=0', 'Non ZLB', 'ZLB')

subplot(2,2,3)
hold on
plot(1:1:hor, zz, 'k-')
plot(1:1:hor, nonzlbirf(:,4), 'b--','LineWidth', 1.5)
plot(1:1:hor, zlbirf(:,4), 'r-o', 'LineWidth', 1.5);
xlabel('quarter')
title('ZLB-dependent: GDP - RZ')
axis tight
legend('y=0', 'Non ZLB', 'ZLB')

subplot(2,2,4)
hold on
plot(1:1:hor, zz, 'k-')
plot(1:1:hor, irfa(2,:), 'b--','LineWidth', 1.5)
plot(1:1:hor, irfb(2,:), 'r-o', 'LineWidth', 1.5);
xlabel('quarter')
title('ZLB-dependent: GDP - my TVAR')
axis tight
legend('y=0', 'Non ZLB', 'ZLB')

% conclusion - we want our TVAR IRFs to be hump-shaped and peak at around
% 0.2-0.3 range

