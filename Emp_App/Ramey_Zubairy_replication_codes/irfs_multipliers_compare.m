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

lin_results = struct(); % struct to store linear model results

% compute the IRFs now with various methods
irfTypes = {'LP', 'LP_BC', 'LP_Penalised', 'LP_Lagaug', 'SVAR', 'VARLP_Avg'};

% standard LP, code from RZ2018
[lin_results.IRF.LP, lin_results.CI.LP]=linlp(data,x,hor,rpos,transformation, clevel, opt, 1); 
liny = lin_results.IRF.LP;
confidencey = lin_results.CI.LP;
% bias correction
[lin_results.IRF.LP_BC, lin_results.CI.LP_BC]=linlp_biascorrect(data,x,hor,rpos,transformation, clevel, opt);
% penalised LP
[lin_results.IRF.LP_Penalised, lin_results.CI.LP_Penalised]=linlp_penalised(data,x,hor,rpos,transformation, clevel, opt, nlag); 
% lag-augmented LP
[lin_results.IRF.LP_Lagaug, lin_results.CI.LP_Lagaug]=linlp_lagaug(data(2:end, :),x_la,hor,rpos,transformation, clevel, opt, 1);
% VAR
[lin_results.IRF.SVAR, lin_results.CI.SVAR] = SVAR(data, x, hor, rpos, transformation, clevel, opt, nlag, 0);
% VAR / LP Averaged (TBD)
[lin_results.IRF.VARLP_Avg, lin_results.CI.VARLP_Avg] = VARLP_Avg(data, x, hor, rpos, transformation, clevel, opt, nlag, 0, 0);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %NON-LINEAR
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear x;
x=xorig;
fu=fu(nlag+1:end);

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

[stateay, stateby, confidenceya, confidenceyb]=statelp_rz(data,x,hor,rpost,transformation, clevel, opt); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run('setup.m') % load settings, colours etc. from this script

% | Figure 5 - the IRFs to GDP and Gov. Spending, Linear Model ONLY |
% -- Plot Gov Spending Response -- %
i=1;
figure(5)
zz = zeros(1, hor);
n = length(irfTypes);
subplot(2,2,1)
plot(1:1:hor, zz, 'k-', 'HandleVisibility', 'off') % plot line at y=0 to show x-axis
hold on

h = zeros(1, n);
for j = 1:length(irfTypes)
    irfType = irfTypes{j};
    h(j) = plot(1:1:hor, lin_results.IRF.(irfType)(i,:), plotStyles{j}{:}, 'DisplayName', irfType);
    hold on;
end

axis tight
ylabel('Government Spending')
lgd = legend('Location', 'northwest');
lgd.Interpreter = 'none'; 

% BELOW NEEDS TO BE FIXED ONCE S.E.s are done

subplot(2,2,2)
grpyat=[(1:1:hor)', confidencey(1,:,i)'; (hor:-1:1)' confidencey(2,hor:-1:1,i)'];
patch(grpyat(:,1), grpyat(:,2), [0.7 0.7 0.7],'edgecolor', [0.7 0.7 0.7]);
hold on
plot(1:1:hor, zz, 'k-', 'HandleVisibility', 'off')
hold on 
plot(1:1:hor, liny(i,:), 'k','LineWidth', 1.5)
title('Linear')
axis tight

% -- Plot GDP Response -- %
i=2;
subplot(2,2,3)
plot(1:1:hor, zz, 'k-', 'HandleVisibility', 'off') % plot line at y=0 to show x-axis
hold on

h = zeros(1, length(irfTypes));
for j = 1:length(irfTypes)
    irfType = irfTypes{j};
    h(j) = plot(1:1:hor, lin_results.IRF.(irfType)(i,:), plotStyles{j}{:}, 'DisplayName', irfType);
    hold on;
end

axis tight
ylabel('GDP')
lgd = legend('Location', 'northwest');
lgd.Interpreter = 'none'; 

% BELOW NEEDS TO BE FIXED ONCE S.E.s are done
subplot(2,2,4)
grpyat=[(1:1:hor)', confidencey(1,:,i)'; (hor:-1:1)' confidencey(2,hor:-1:1,i)'];
patch(grpyat(:,1), grpyat(:,2), [0.7 0.7 0.7],'edgecolor', [0.7 0.7 0.7]);
hold on
plot(1:1:hor, zz, 'k-', 'HandleVisibility', 'off')
hold on 
plot(1:1:hor, liny(i,:), 'k','LineWidth', 1.5)
title('Linear')
axis tight

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Multipliers: HAVE EDITED LINEAR ONLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute linear multipliers
for j = 1:length(irfTypes)
    irfType = irfTypes{j};
    lin_results.cum_mult.(irfType) = cumsum(lin_results.IRF.(irfType)(2, :)) ./ cumsum(lin_results.IRF.(irfType)(1, :));
end

cum_mult_statea= cumsum(stateay(2,:))./cumsum(stateay(1,:));
cum_mult_stateb= cumsum(stateby(2,:))./cumsum(stateby(1,:));

% reads the multiplier standard errors that are computed via STATA
if shockchoice==1 % newsy
    std1=xlsread('Multiplier-Standard-Errors.xlsx',2);
    stdlin=std1(1:20,3);
    if statechoice==1 % slack
        std2=xlsread('Multiplier-Standard-Errors.xlsx',2);
        stdb=std2(1:20,5);
        stda=std2(1:20,7);
    elseif statechoice==2 % ZLB
        std2=xlsread('Multiplier-Standard-Errors.xlsx',3);
        stdb=std2(1:20,5);
        stda=std2(1:20,7);
    end
else
    std1=xlsread('Multiplier-Standard-Errors.xlsx',4);
    stdlin=std1(1:20,3);
    if statechoice==1
        std2=xlsread('Multiplier-Standard-Errors.xlsx',4);
        stdb=std2(1:20,5);
        stda=std2(1:20,7);
    elseif statechoice==2
        std2=xlsread('Multiplier-Standard-Errors.xlsx',5);
        stdb=std2(1:20,5);
        stda=std2(1:20,7);
    end
end

multconfb=[cum_mult_stateb+clevel*stdb'; cum_mult_stateb-clevel*stdb'];
multconfa=[cum_mult_statea+clevel*stda'; cum_mult_statea-clevel*stda'];

% the linear model multiplier confidence intervals are placeholders 
for j = 1:length(irfTypes)
    irfType = irfTypes{j};
    lin_results.cum_mult_CI.(irfType) = [lin_results.cum_mult.(irfType)+clevel*stdlin'; lin_results.cum_mult.(irfType)-clevel*stdlin'];
end

% Figure for multipliers
figure(6)

for i = 1:length(irfTypes)
    irfType = irfTypes{i};
    subplot(2,length(irfTypes),i)
    
    grpyat=[(1:1:hor)', lin_results.cum_mult_CI.(irfType)(1,:)'; (hor:-1:1)', lin_results.cum_mult_CI.(irfType)(2,hor:-1:1)']; % create matrix for the CIs (basically it ends up tracing 1 -> hor upper bounds, then back around hor -> 1 lower bounds)
    patch(grpyat(:,1), grpyat(:,2), patch_colors{i},'edgecolor', patch_colors{i}, 'DisplayName', strcat(irfType, ' CI')); % plot gray patch for CIs
    hold on 
    plot(1:1:hor, lin_results.cum_mult.(irfType), 'Color', line_colors{i}, 'LineWidth', 1.5, 'DisplayName', irfType); % plot point estimates for multiplier

    axis tight
    title('Linear: cumulative spending multiplier');
    legend('Location', 'northeast')
end

% state dependent plot
subplot(2,length(irfTypes),5)
grpyat=[(1:1:hor)', multconfa(1,:)'; (hor:-1:1)' multconfa(2,hor:-1:1)'];
patch(grpyat(:,1), grpyat(:,2), [0.7 0.7 0.7],'edgecolor', [0.7 0.7 0.7]);
hold on 
plot(1:1:hor, cum_mult_statea, 'b--', 'LineWidth', 1.5);
hold on
plot(1:1:hor, cum_mult_stateb+clevel*stdb', 'r--', 1:1:hor, cum_mult_stateb-clevel*stdb', 'r--','LineWidth', 1)
hold on
plot(1:1:hor, cum_mult_stateb, 'r-o', 'LineWidth', 1.5);
title('State dependent: cumulative spending multiplier')
xlabel('quarter')
axis tight



