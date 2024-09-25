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

%% analytical var results from STATA
stata_irfs = readtable("01_Test_NOTFORUSE/stata_lin_svar.csv")

%% test linear IRFs
% MATLAB irf function with MC standard errors / CIs
horizon = 21;
y = [x(:,rpos), data];
p = 4;
[T, k] = size(y);

model = varm(k,p);
[mdl, ~, ~, res] = estimate(model, y);
[irf_var, lower, upper] = irf(mdl, E=res, NumObs = horizon, NumPaths=200); 

% STATA 
shockind = 1;
responseind = 2;

%% linear IRFs with bootstrapping (manual code)
bootstrap = 1;
B = 200;
rind = [2 3];
sind = 1;
normalise = 0;
[myirf_var, myirfvar_ci] = linSVAR(y,horizon,nlag,rind,sind,clevel,bootstrap,B,normalise);

%%
close all
fig = figure;
ax1 = subplot(2,3,1);
hold on
plot(irf_var(:,shockind,2), 'LineStyle','-')
plot(lower(:,shockind,2), 'LineStyle','--')
plot(upper(:,shockind,2), 'LineStyle','--')
legend('irf', 'lower', 'upper')
ylabel('IRF to Gov Spend')
title('MATLAB IRF SEs using Monte Carlo CIs')

ax2 = subplot(2,3,2);
hold on
plot(stata_irfs.g_oirf, 'LineStyle','-')
plot(stata_irfs.g_lower, 'LineStyle','--')
plot(stata_irfs.g_upper, 'LineStyle','--')
legend('irf', 'lower', 'upper')
ylabel('IRF to Gov Spend')
title('STATA (RZ method) IRF SEs using Analytical CIs')
axis tight

ax3 = subplot(2,3,3);
hold on
plot(myirf_var(1,:), 'LineStyle','-')
plot(squeeze(myirfvar_ci(1,1,:)), 'LineStyle','--')
plot(squeeze(myirfvar_ci(2,1,:)), 'LineStyle','--')
legend('irf', 'lower', 'upper')
ylabel('IRF to Gov Spend')
title('MATLAB Custom Code for Lutkepohl (2000) Bootstrap Method')
axis tight

ax4 = subplot(2,3,4);
hold on
plot(irf_var(:,shockind,3), 'LineStyle','-')
plot(lower(:,shockind,3), 'LineStyle','--')
plot(upper(:,shockind,3), 'LineStyle','--')
legend('irf', 'lower', 'upper')
ylabel('IRF to GDP')

ax5 = subplot(2,3,5);
hold on
plot(stata_irfs.y_oirf, 'LineStyle','-')
plot(stata_irfs.y_lower, 'LineStyle','--')
plot(stata_irfs.y_upper, 'LineStyle','--')
legend('irf', 'lower', 'upper')
ylabel('IRF to GDP')
axis tight

ax6 = subplot(2,3,6);
hold on
plot(myirf_var(2,:), 'LineStyle','-')
plot(squeeze(myirfvar_ci(1,2,:)), 'LineStyle','--')
plot(squeeze(myirfvar_ci(2,2,:)), 'LineStyle','--')
legend('irf', 'lower', 'upper')
ylabel('IRF to GDP')
axis tight

linkaxes([ax1, ax2, ax3, ax4, ax5, ax6], 'y')
linkaxes([ax1, ax2, ax3, ax4, ax5, ax6], 'x')

fig.Position = [100 100 1200 900];  % [left bottom width height]

% Save the figure with high resolution
print(fig, 'Fig/linearSVAR_mc_bs_analytical', '-dpng', '-r300');

%% plot difference in standard errors
stata_seg = abs(stata_irfs.g_upper - stata_irfs.g_oirf) / clevel;
stata_sey = abs(stata_irfs.y_upper - stata_irfs.y_oirf) / clevel;

bs_seg = abs(myirf_var(1,:) - squeeze(myirfvar_ci(2,1,:))')/clevel;
bs_sey = abs(myirf_var(2,:) - squeeze(myirfvar_ci(2,2,:))')/clevel;

% standard error discrepancy
[(stata_seg - bs_seg'), (stata_sey - bs_sey')]
[(stata_seg - bs_seg')./stata_irfs.g_oirf, (stata_sey - bs_sey')./stata_irfs.y_oirf]

%%

% lpw var script
[irf_var_lpw, ~] = linSVAR(data, x, hor, rpos, transformation, clevel, opt, nlag, 0);

% lp comparison
[irf_lp, ~] = linlp(data,x,hor,rpos,transformation, clevel, opt, nlag, 0);

% plot
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
sgtitle('Linear Model IRF to Gov Spend - VAR MATLAB, Manual VAR - LPW, LP')

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


%% FIRST, TEST OUR TVAR and corr SEs (SHOWED THIS CODE WORKS)
% Note: in RZ TVAR figs, the blue represents non-recession / non-zlb, and
% red means reces
run('figure_irfs_multipliers_setup.m')

y = y_tvar; % ordered as shock, gov, gdp
[t,k] = size(y);
I = fu;
p=4;
hor=20;
rind=[2,3]; % irfs to gov and GDP
sind=1; % shock ordered first in the y vector
B = 100;

[irfa, irfb, ~, ~, ~, ~, beta] = stateSVAR(y, I, p, hor,rind,sind, B);

% check eig;
beta_a = beta(:,:,1);
beta_b = beta(:,:,2);

ones_comp = [eye(k*p) zeros(k*p,k)];
ones_comp = ones_comp(1:end-k, 1:end-k);
companion_matrix_a = [beta_a(2:end,:)'; ones_comp];
companion_matrix_b = [beta_b(2:end,:)'; ones_comp];


ind = [1];
for i=1:k
    k_temp = [];
    for j = 1:p
        k_temp = [k_temp i+j+1];
    end
    ind = [ind k_temp];
end

if any(abs(eig(companion_matrix_a)) >= 1)
    disp('A may be explosive')
end
if any(abs(eig(companion_matrix_b)) >= 1)
    disp('B may be explosive')
end

figure;
subplot(1,2,1);
plot(irfa(1,:))
hold on
plot(irfb(1,:))
legend('State A', 'State B - Non')
title('Shock to Gov Spend')

subplot(1,2,2);
plot(irfa(2,:))
hold on
plot(irfb(2,:))
legend('State A', 'State B - Non')
title('Shock to GDP')

%% State dependent LP
run('figure_irfs_multipliers_setup.m')

x=xorig;
In=fu(nlag+1:end); % truncated indicator variable for LP regression

if trend==4 % allow quartic trend
    x=[t', tsq', tcu', tqu', constant, (1-In).*constant, In.*shock, (1-In).*shock, repmat(In,1,size(x,2)).*x, repmat((1-In),1,size(x,2)).*x];
    rpost=7; 
elseif trend==2 % allow quadratic trend
    x=[t', tsq', constant, (1-In).*constant, In.*shock, (1-In).*shock, repmat(In,1,size(x,2)).*x, repmat((1-In),1,size(x,2)).*x];
    rpost=5; 
elseif trend==0 %no trend
    x=[constant, (1-In).*constant, In.*shock, (1-In).*shock, repmat(In,1,size(x,2)).*x, repmat((1-In),1,size(x,2)).*x];
    rpost=3; % position of shock
end

% settings
bootstrap = 1;
nstraps=99;
emp=0;

[stateay, stateby, confidenceya, ~, ~, ~] = statelp(data,x,hor,rpost,transformation, clevel, opt, bootstrap, emp, p, nstraps);

disp(confidenceya(1,:,1)')
%% model averaged estimator - point

lambda = 0.5;
% first do it manually
avgirfa = lambda * irfa + (1-lambda) * stateay;
avgirfb = lambda * irfb + (1-lambda) * stateby;

% settings for stateSVARLP_avg
B = 100;
% no bootstrap argument as only possibility is bootstrap

[mdlavga, mdlavgb, ~, ~, ~, ~] = stateSVARLP_avg(hor, nlag, B, lambda, ...
    data,x,rpost,transformation, clevel, opt, ...
    y, I, rind, sind);

%% Compare the three - validated :)
d=1;
figure;
subplot(1,4,1)
hold on
plot(0:19, irfa(d,:), 'Color', 'red')
plot(0:19, irfb(d,:), 'Color', 'blue')
legend('A - VAR', 'B - VAR') 

subplot(1,4,2)
hold on
plot(0:19, stateay(d,:), 'Color', 'red', 'LineStyle','--')
plot(0:19, stateby(d,:), 'Color', 'blue', 'LineStyle','--')
legend('A - LP', 'B - LP')

subplot(1,4,3)
hold on
plot(0:19, avgirfa(d,:), 'Color', 'red', 'LineStyle','--')
plot(0:19, avgirfb(d,:), 'Color', 'blue', 'LineStyle','--')
legend('A - AvgMan', 'B - AvgMan')

subplot(1,4,4)
hold on
plot(0:19, mdlavga(d,:), 'Color', 'red', 'LineStyle','--')
plot(0:19, mdlavgb(d,:), 'Color', 'blue', 'LineStyle','--')
legend('A - AvgFunc', 'B - AvgFunc')


%% model averaged estimator - bootstrap / p val
lambda = 0.5;
B = 500;

[mdlavga, mdlavgb, ~, ~, bs_beta_dist, pval] = stateSVARLP_avg(hor, nlag, B, lambda, ...
    data,x,rpost,transformation, clevel, opt, ...
    y, I, rind, sind);

pvala = squeeze(pval(1,:,:));
pvalb = squeeze(pval(2,:,:));

a = array2table(pvala', 'VariableNames',{'Gov Spending', 'GDP'});
b = array2table(pvalb', 'VariableNames',{'Gov Spending', 'GDP'});
%%
output = evalc('disp(array2table(table2array(a), ''VariableNames'', a.Properties.VariableNames))');
disp(output);

%%
bs_beta_dista = squeeze(bs_beta_dist(1,:,:,:));
bs_beta_distb = squeeze(bs_beta_dist(2,:,:,:));
ri = 1;
hi = 4;

close all
figure;
hist(bs_beta_dista(:,ri,hi), 20)
hold on
vline(mdlavga(ri, hi))




% conclusion - we want our TVAR IRFs to be hump-shaped and peak at around
% 0.2-0.3 range

