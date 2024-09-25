%%% SCRIPT DESIGNED TO TEST STATE DEPENDENT VAR BOOTSTRAP

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

%% test TVAR vs RZ [crucially, both point and CIs]
% FIRST, TEST OUR TVAR (SHOWED THIS CODE WORKS)
y = y_tvar; % ordered as shock, gov, gdp
I = fu;
p=4;
hor=20;
B = 1000;
rind=[2,3]; % irfs to gov and GDP
sind=1; % shock ordered first in the y vector
bootstrap = 1;
method = "normal";

[irfa, irfb, cia, cib, ~, ~, beta] = stateSVAR(y, I, p, hor,rind,sind, B, bootstrap, clevel, method);

%%
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

fig = figure;

ri = 1;
ax1 = subplot(2,2,1);
hold on
plot(1:1:hor, zz, 'k-', 'HandleVisibility','off')
plot(1:1:hor, nonzlbirf(:,1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Non ZLB: Point')
plot(1:1:hor, nonzlbirf(:,2), 'b--', 'LineWidth', 1, 'DisplayName', 'Non ZLB: CI')
plot(1:1:hor, nonzlbirf(:,3), 'b--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot(1:1:hor, zlbirf(:,1), 'r-', 'LineWidth', 1.5, 'DisplayName', 'ZLB: Point')
plot(1:1:hor, zlbirf(:,2), 'r--', 'LineWidth', 1, 'DisplayName', 'ZLB: CI')
plot(1:1:hor, zlbirf(:,3), 'r--', 'LineWidth', 1, 'HandleVisibility', 'off')
xlabel('quarter')
title('ZLB-dependent: Government spending - RZ')
axis tight
legend('show')

ax2 = subplot(2,2,2);
hold on
plot(1:1:hor, zz, 'k-', 'HandleVisibility','off')
plot(1:1:hor, irfa(ri,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Non ZLB: Point')
plot(1:1:hor, squeeze(cia(1,ri,:)), 'b--', 'LineWidth', 1, 'DisplayName', 'Non ZLB: CI')
plot(1:1:hor, squeeze(cia(2,ri,:)), 'b--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot(1:1:hor, irfb(ri,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'ZLB: Point')
plot(1:1:hor, squeeze(cib(1,ri,:)), 'r--', 'LineWidth', 1, 'DisplayName', 'ZLB: CI')
plot(1:1:hor, squeeze(cib(2,ri,:)), 'r--', 'LineWidth', 1, 'HandleVisibility', 'off')
xlabel('quarter')
title('ZLB-dependent: Government spending - self coded bootstrap TVAR')
axis tight
legend('show')

ri = 2;
ax3 = subplot(2,2,3);
hold on
plot(1:1:hor, zz, 'k-', 'HandleVisibility','off')
plot(1:1:hor, nonzlbirf(:,4), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Non ZLB: Point')
plot(1:1:hor, nonzlbirf(:,5), 'b--', 'LineWidth', 1, 'DisplayName', 'Non ZLB: CI')
plot(1:1:hor, nonzlbirf(:,6), 'b--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot(1:1:hor, zlbirf(:,4), 'r-', 'LineWidth', 1.5, 'DisplayName', 'ZLB: Point')
plot(1:1:hor, zlbirf(:,5), 'r--', 'LineWidth', 1, 'DisplayName', 'ZLB: CI')
plot(1:1:hor, zlbirf(:,6), 'r--', 'LineWidth', 1, 'HandleVisibility', 'off')
xlabel('quarter')
title('ZLB-dependent: GDP - RZ')
axis tight
legend('show')

ax4 = subplot(2,2,4);
hold on
plot(1:1:hor, zz, 'k-', 'HandleVisibility','off')
plot(1:1:hor, irfa(ri,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Non ZLB: Point')
plot(1:1:hor, squeeze(cia(1,ri,:)), 'b--', 'LineWidth', 1, 'DisplayName', 'Non ZLB: CI')
plot(1:1:hor, squeeze(cia(2,ri,:)), 'b--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot(1:1:hor, irfb(ri,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'ZLB: Point')
plot(1:1:hor, squeeze(cib(1,ri,:)), 'r--', 'LineWidth', 1, 'DisplayName', 'ZLB: CI')
plot(1:1:hor, squeeze(cib(2,ri,:)), 'r--', 'LineWidth', 1, 'HandleVisibility', 'off')
xlabel('quarter')
title('ZLB-dependent: GDP - self coded bootstrap TVAR')
axis tight
legend('show')

linkaxes([ax1, ax2, ax3, ax4], 'y')

fig.Position = [100 100 800 600];  % [left bottom width height]

saveas(fig, 'Fig/stateSVAR_bs_analytical.png')