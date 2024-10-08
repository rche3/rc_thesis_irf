%%% SCRIPT DESIGNED TO TEST STATE DEPENDENT LP BOOTSTRAP

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

%% test SD LP bootstrap vs RZ [crucially, both point and CIs], nlag
% FIRST, TEST OUR TVAR (SHOWED THIS CODE WORKS)
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

bootstrap = 1;
nstraps = 100;
method = 'normal'
% Second, use DW bootstrap
[myirfa, myirfb, mycia, mycib, ~,~] = statelp(data,shock,z,In,hor, ...
    transformation,clevel,opt, ...
    bootstrap,method,nlag,nstraps);

%%
close all
shockchoice=1; % 1 means news shock; 2 means BP shock

hor=20;

fig = figure;
zz=zeros(1,hor);

ri = 1;
ax1 = subplot(2,2,1);
hold on
plot(1:1:hor, zz, 'k-', 'HandleVisibility','off')
plot(1:1:hor, rzirfa(ri,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'ZLB: Point')
plot(1:1:hor, squeeze(rzcia(1,ri,:)), 'b--', 'LineWidth', 1, 'DisplayName', 'ZLB: CI')
plot(1:1:hor, squeeze(rzcia(2,ri,:)), 'b--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot(1:1:hor, rzirfb(ri,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Non ZLB: Point')
plot(1:1:hor, squeeze(rzcib(1,ri,:)), 'r--', 'LineWidth', 1, 'DisplayName', 'Non ZLB: CI')
plot(1:1:hor, squeeze(rzcib(2,ri,:)), 'r--', 'LineWidth', 1, 'HandleVisibility', 'off')
ylabel('Gov Spending Response')
xlabel('quarter')
title('ZLB-dependent: Government spending - RZ Original')
axis tight
legend('show')

ax2 = subplot(2,2,2);
hold on
plot(1:1:hor, zz, 'k-', 'HandleVisibility','off')
plot(1:1:hor, myirfa(ri,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'ZLB: Point')
plot(1:1:hor, squeeze(mycia(1,ri,:)), 'b--', 'LineWidth', 1, 'DisplayName', 'ZLB: CI')
plot(1:1:hor, squeeze(mycia(2,ri,:)), 'b--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot(1:1:hor, myirfb(ri,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Non ZLB: Point')
plot(1:1:hor, squeeze(mycib(1,ri,:)), 'r--', 'LineWidth', 1, 'DisplayName', 'Non ZLB: CI')
plot(1:1:hor, squeeze(mycib(2,ri,:)), 'r--', 'LineWidth', 1, 'HandleVisibility', 'off')
ylabel('Gov Spending Response')
xlabel('quarter')
title('ZLB-dependent: Government spending - Self coded Bootstrap')
axis tight
legend('show')

% GDP
ri = 2;

ax3 = subplot(2,2,3);
hold on
plot(1:1:hor, zz, 'k-', 'HandleVisibility','off')
plot(1:1:hor, rzirfa(ri,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'ZLB: Point')
plot(1:1:hor, squeeze(rzcia(1,ri,:)), 'b--', 'LineWidth', 1, 'DisplayName', 'ZLB: CI')
plot(1:1:hor, squeeze(rzcia(2,ri,:)), 'b--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot(1:1:hor, rzirfb(ri,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Non ZLB: Point')
plot(1:1:hor, squeeze(rzcib(1,ri,:)), 'r--', 'LineWidth', 1, 'DisplayName', 'Non ZLB: CI')
plot(1:1:hor, squeeze(rzcib(2,ri,:)), 'r--', 'LineWidth', 1, 'HandleVisibility', 'off')
ylabel('GDP Response')
xlabel('quarter')
title('ZLB-dependent: GDP - RZ Original')
axis tight
legend('show')

ax4 = subplot(2,2,4);
hold on
plot(1:1:hor, zz, 'k-', 'HandleVisibility','off')
plot(1:1:hor, myirfa(ri,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'ZLB: Point')
plot(1:1:hor, squeeze(mycia(1,ri,:)), 'b--', 'LineWidth', 1, 'DisplayName', 'ZLB: CI')
plot(1:1:hor, squeeze(mycia(2,ri,:)), 'b--', 'LineWidth', 1, 'HandleVisibility', 'off')
plot(1:1:hor, myirfb(ri,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Non ZLB: Point')
plot(1:1:hor, squeeze(mycib(1,ri,:)), 'r--', 'LineWidth', 1, 'DisplayName', 'Non ZLB: CI')
plot(1:1:hor, squeeze(mycib(2,ri,:)), 'r--', 'LineWidth', 1, 'HandleVisibility', 'off')

ylabel('GDP Response')
xlabel('quarter')
title('ZLB-dependent: GDP - Self Coded Bootstrap')
axis tight
legend('show')

linkaxes([ax1, ax2], 'y')
linkaxes([ax3, ax4], 'y')
fig.Position = [100 100 800 600];  % [left bottom width height]

saveas(fig, 'Fig/statelp_bs_analytical.png')