clear all
close all
%clc
warning off;

addpath(genpath('Plots'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT SPECIFICATION CHOICES / SETUP
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%NON-LINEAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nlr = struct(); % struct to store estimates for non linear results (nlr)

% compute the IRFs now with various methods
% irfTypes = {'LP', 'LP_BC', 'LP_Penalised', 'LP_Lagaug', 'SVAR', 'VARLP_Avg'};
irfTypes = {'LP', 'SVAR', 'SVARLPAvg'};

L = length(irfTypes);

%% Compute IRFs
% key settings
emp = 0;
p = nlag;
lambda = 0.5;
nstraps = 200;

e=1; est=irfTypes{e}; 
% STATE DEPENDENT LP
x=xorig;
In=fu(nlag+1:end);
bootstrap = 0;
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
[nlr.IRFa.(est), nlr.IRFb.(est), nlr.CIa.(est), nlr.CIb.(est),~, nlr.pval.(est)]=statelp(data,x,hor,rpost,transformation, clevel, opt, bootstrap, emp, p, nstraps); 

e=2; est = irfTypes{e};
% SVAR
I = fu;
y = y_tvar;
p = nlag;
rind = [2:3];
sind = 1;
bootstrap = 1;
[nlr.IRFa.(est), nlr.IRFb.(est),nlr.CIa.(est), nlr.CIb.(est), ~, nlr.pval.(est), ~] = stateSVAR(y,I,p,hor,rind,sind, nstraps, bootstrap, clevel);

e=3; est = irfTypes{e};
% SVAR / LP MODEL AVG
[nlr.IRFa.(est), nlr.IRFb.(est), nlr.CIa.(est), nlr.CIb.(est), nlr.bsdist.(est), nlr.pval.(est)] = stateSVARLP_avg(hor, nlag, nstraps, lambda, ...
    data,x,rpost,transformation, clevel, opt, ...
    y, I, rind, sind);

%% Figure 5 - Non Linear IRFs
close all
zz=zeros(1,hor);
f5 = figure;
ylabels = {'Government Spending', 'GDP'};

for j=1:L
    est = irfTypes{j};
    for i=1:2
        subplot(2,L,j+L*(i-1))
        plot(1:1:hor, zz, 'k-', 'HandleVisibility','off')
        hold on
        plot(1:1:hor, nlr.IRFa.(est)(i,:), 'b-', ...
             1:1:hor, squeeze(nlr.CIa.(est)(1,i,:)), 'b--', ...
             1:1:hor, squeeze(nlr.CIa.(est)(2,i,:)), 'b--');
        plot(1:1:hor, nlr.IRFb.(est)(i,:), 'r-', ...
             1:1:hor, squeeze(nlr.CIb.(est)(1,i,:)), 'r--', ...
             1:1:hor, squeeze(nlr.CIb.(est)(2,i,:)), 'r--');
        axis tight
        ylabel(ylabels{i})
        xlabel('Horizon (Quarters)')
        title(est)
        legend('A: Point', 'A: 95% CI Lower', 'A: 95% CI Upper', 'B: Point', 'B: 95% CI Lower', 'B: 95% CI Upper')
    end
end
sgtitle('State Dependent IRF Estimates')

set(f5, 'Position', [100, 100, 1200, 800])
% saveas(f5,'fig/lp_var_mdlavg_point_est.png')

%% P value table
state = 1; % 1 for State A, 2 for State B

est = 'SVARLPAvg';

pval = squeeze(nlr.pval.(est)(state,:,:));
a = array2table(pval', 'VariableNames',{'Gov Spending', 'GDP'})

fig = figure;
subplot(1,2,1)
bar(pval(1,:))
legend('Gov Spend')
subplot(1,2,2)
bar(pval(2,:))
legend('GDP')
sgtitle(strcat('Two Tailed P Values vs Horizon for:', est));

% saveas(fig, 'fig/SVARLP_Pval_BS.png')

%% Look at distribution of bootstrap irfs for SVAR LP Avg
est = 'SVARLPAvg';
state = 1; % 1 for State A, 2 for State B
temphor = 1;
tempresp = 1;

bs_irf_dist = nlr.bsdist.(est)(state,:,tempresp, temphor);
point_est = nlr.IRFa.SVARLPAvg(tempresp, temphor);

figure;
hist(bs_irf_dist, 20)
hold on
vline(point_est)



