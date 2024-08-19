%% Setup
clear all
close all
%clc
warning off;

cd /Users/rogerchen/Documents/MATLAB/rc_thesis_irf/Emp_App/RZ_Active_Rep/

%% Input Specs
timesample=1; % 1 means 1889-2015; 2 means post WWII
statechoice=1;  % 1 means default unemp threshold; 2 means default ZLB threshold
shockchoice=1; % 1 means news shock; 2 means BP shock
transformation=1; % 1 means Gordon-Krenn; 2 means Hall-Barro-Redlick
taxchoice=0; %0 means no taxes as control; 1 means ad taxy as a control
datafirstdiff=0; % means to difference the key GDP / Gov Spending series 

trend=0; %4 means quartic; 2 means quadratic; 0 means no trend  
nlag=4;% number of lags of control variables
opt=0; % 0 means Newey-West with bandwidth=i; 1 means optimal bandwidth (takes much longer to run)'

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

hor=20; % horizon for IRFs
clevel= 1.96; % confidence level % 1.96  1.64;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run('figure_irfs_multipliers_setup.m')

%% Estimating Linear Model w/ Analytical SEs
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

% temporarily redefine data
% T_ss = 2500;
% data = randn(T_ss, 2);
% x = [ones(T_ss, 1) randn(T_ss, 13)];

% lag_regressor = lagmatrix(regressor, [1:4]); 
% data = data(5:end, :);
% x = [ones(T_ss-4, 1) regressor(5:end) lag_regressor(5:end, :)];

% Generate the standard LP IRFs
[liny, confidencey]=linlp_rc(data,x,hor,rpos,transformation,clevel,opt,0,nlag);
se = abs(liny(1, :) - confidencey(1, :, 1))/clevel;
se = se';

%% Compute Bootstrap SEs
clear bs_confidencey
[liny, bs_confidencey, bs_beta_means] = linlp_rc(data,x,hor,rpos,transformation, clevel, opt, 1, nlag, 0);
bs_se = abs(liny(1, :) - bs_confidencey(1, :, 1))/clevel;
bs_se = bs_se';

% %% Override Bootstrap SEs with Lag-Augmented CIs
% clear bs_confindencey
% [liny, bs_confidencey, bs_beta_means] = linlp_rc(data(2:end, :),x_la,hor,rpos,transformation, clevel, opt, 1, nlag+1);
% bs_se = abs(liny(1, :) - bs_confidencey(1, :, 1))/clevel';
% bs_se = bs_se';
% 
% %% Override Bootstrap SEs with sample mean of NWest Corrected
% clear bs_confindencey
% [liny, bs_confidencey, bs_beta_means] = linlp_rc(data,x,hor,rpos,transformation, clevel, opt, 1, nlag, 1);
% bs_se = abs(liny(1, :) - bs_confidencey(1, :, 1))/clevel';
% bs_se =on bs_se';

% Plot Bootstrap SEs and Analytical SEs
close all
zz=zeros(1,hor);
i=1;
figure(1)
subplot(2,1,1)
hold on
% plot CIs - analytical
grpyat=[(1:1:hor)', confidencey(1,:,i)'; (hor:-1:1)' confidencey(2,hor:-1:1,i)'];
patch(grpyat(:,1), grpyat(:,2), [0.5 0.5 0.5],'edgecolor', [0.5 0.5 0.5], 'FaceAlpha', 0.5);
% plot CIs - bootstrapped
grpyat=[(1:1:hor)', bs_confidencey(1,:,i)'; (hor:-1:1)' bs_confidencey(2,hor:-1:1,i)'];
patch(grpyat(:,1), grpyat(:,2), [0.5 0.9 0.5],'edgecolor', [0.5 0.9 0.5], 'FaceAlpha', 0.3);
% plot point estimates
plot(1:1:hor, liny(i,:), 'Color', [0.08 0.64 0.4],'LineWidth', 1.5)
plot(1:1:hor, zz, 'k-')
axis tight
ylabel('Government spending')
xlabel('quarter')
legend('Analytical CI', 'BS CI', 'Point Estimates')

i=2;
subplot(2,1,2)
hold on
% plot CIs - analytical
grpyat=[(1:1:hor)', confidencey(1,:,i)'; (hor:-1:1)' confidencey(2,hor:-1:1,i)'];
patch(grpyat(:,1), grpyat(:,2), [0.5 0.5 0.5],'edgecolor', [0.5 0.5 0.5], 'FaceAlpha', 0.5);
% plot CIs - bootstrapped
grpyat=[(1:1:hor)', bs_confidencey(1,:,i)'; (hor:-1:1)' bs_confidencey(2,hor:-1:1,i)'];
patch(grpyat(:,1), grpyat(:,2), [0.5 0.5 0.9],'edgecolor', [0.5 0.5 0.9], 'FaceAlpha', 0.3);
% plot ploint estimates
plot(1:1:hor, liny(i,:), 'Color', 'black','LineWidth', 1.5)
plot(1:1:hor, zz, 'k-')
axis tight
ylabel('GDP')
xlabel('quarter')
legend('Analytical CI', 'BS CI', 'Point Estimates')


%% Check convergence of bootstrap means

i=1;
figure(3)
subplot(2,1,1)
hold on
% plot point estimates
plot(1:1:hor, liny(i,:), 'Color', 'red','LineWidth', 1.5)
plot(1:1:hor, bs_beta_means(i,:), 'Color', 'blue','LineWidth', 1.5)
plot(1:1:hor, zz, 'k-')
axis tight
ylabel('Government spending')
legend('LP IRF', 'BS IRF')

i=2;
subplot(2,1,2)
hold on
plot(1:1:hor, liny(i,:), 'Color', 'red','LineWidth', 1.5)
plot(1:1:hor, bs_beta_means(i,:), 'Color', 'blue','LineWidth', 1.5)
plot(1:1:hor, zz, 'k-')
axis tight
ylabel('GDP')
xlabel('quarter')
legend('LP IRF', 'BS IRF')


%% NON-LINEAR Model
clear x;
x=xorig;
fu=fu(nlag+1:end); % fu = 1 if the previous period exceeded threshold

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

[stateay, stateby, confidenceya, confidenceyb]=statelp_rz(data,x,hor,rpost,transformation, clevel, opt, 0); 

%% NON-LINEAR Model - Bootstrapped CIs

[stateay, stateby, bs_confidenceya, bs_confidenceyb] = statelp_rz(data,x,hor,rpost,transformation, clevel, opt, 0); 
%% Plotting Non-Linear
close all
zz=zeros(1,hor);
i=1;
figure(1)
tiledlayout(2,2);
nexttile;
hold on
% plot high unemployment (state A) point estimates and CI
plot(1:1:hor, stateay(i,:), 'Color', [0.08 0.64 0.4],'LineWidth', 1.5)
grpyat=[(1:1:hor)', confidenceya(1,:,i)'; (hor:-1:1)' confidenceya(2,hor:-1:1,i)'];
patch(grpyat(:,1), grpyat(:,2), [0.5 0.5 0.5],'edgecolor', [0.5 0.5 0.5], 'FaceAlpha', 0.5);
% plot CIs - bootstrapped
% grpyat=[(1:1:hor)', bs_confidencey(1,:,i)'; (hor:-1:1)' bs_confidencey(2,hor:-1:1,i)'];
% patch(grpyat(:,1), grpyat(:,2), [0.5 0.9 0.5],'edgecolor', [0.5 0.9 0.5], 'FaceAlpha', 0.3);
plot(1:1:hor, zz, 'k-')
axis tight
ylabel('Government spending')
xlabel('quarter')
title('High unemployment IRFs')

nexttile;
hold on
% plot high unemployment (state A) point estimates and CI
plot(1:1:hor, stateby(i,:), 'Color', [0.08 0.64 0.4],'LineWidth', 1.5)
grpyat=[(1:1:hor)', confidenceyb(1,:,i)'; (hor:-1:1)' confidenceyb(2,hor:-1:1,i)'];
patch(grpyat(:,1), grpyat(:,2), [0.5 0.5 0.5],'edgecolor', [0.5 0.5 0.5], 'FaceAlpha', 0.5);
% plot CIs - bootstrapped
plot(1:1:hor, zz, 'k-')
axis tight
ylabel('GDP')
xlabel('quarter')
title('Low unemployment IRFs')

i=2;
nexttile;
hold on
% plot high unemployment (state A) point estimates and CI
plot(1:1:hor, stateay(i,:), 'Color', [0.08 0.64 0.4],'LineWidth', 1.5)
grpyat=[(1:1:hor)', confidenceya(1,:,i)'; (hor:-1:1)' confidenceya(2,hor:-1:1,i)'];
patch(grpyat(:,1), grpyat(:,2), [0.5 0.5 0.5],'edgecolor', [0.5 0.5 0.5], 'FaceAlpha', 0.5);
% plot CIs - bootstrapped
% grpyat=[(1:1:hor)', bs_confidencey(1,:,i)'; (hor:-1:1)' bs_confidencey(2,hor:-1:1,i)'];
% patch(grpyat(:,1), grpyat(:,2), [0.5 0.9 0.5],'edgecolor', [0.5 0.9 0.5], 'FaceAlpha', 0.3);
plot(1:1:hor, zz, 'k-')
axis tight
ylabel('Government spending')
xlabel('quarter')
title('High unemployment IRFs')

%%% 
nexttile;
hold on
% plot high unemployment (state A) point estimates and CI
plot(1:1:hor, stateby(i,:), 'Color', [0.08 0.64 0.4],'LineWidth', 1.5)
grpyat=[(1:1:hor)', confidenceyb(1,:,i)'; (hor:-1:1)' confidenceyb(2,hor:-1:1,i)'];
patch(grpyat(:,1), grpyat(:,2), [0.5 0.5 0.5],'edgecolor', [0.5 0.5 0.5], 'FaceAlpha', 0.5);
% plot CIs - bootstrapped
plot(1:1:hor, zz, 'k-')
axis tight
ylabel('GDP')
xlabel('quarter')
title('Low unemployment IRFs')
