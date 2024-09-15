clear all
close all
%clc
warning off;

addpath(genpath('Plots'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT SPECIFICATION CHOICES / SETUP
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

nlr = struct(); % struct to store estimates for non linear results (nlr)

% compute the IRFs now with various methods
% irfTypes = {'LP', 'LP_BC', 'LP_Penalised', 'LP_Lagaug', 'SVAR', 'VARLP_Avg'};
irfTypes = {'LP', 'SVAR'};
L = length(irfTypes);

%% Compute IRFs
bootstrap = 0;
nstraps = 199;
emp = 0;
p = nlag;
e=1; est=irfTypes{e};
[nlr.IRFa.(est), nlr.IRFb.(est), nlr.CIa.(est), nlr.CIb.(est), ~, ~, ~, ~]=statelp(data,x,hor,rpost,transformation, clevel, opt, bootstrap, emp, p, nstraps); 

% prepare data for SVAR / SVAR model avg estimation
% e=2; est = irfTypes{e};
% [nlr.IRFa.(est), nlr.IRFb.(est), ~, ~, nlr.pvala.(est), nlr.pvalb.(est)] = stateSVAR(ya,yb,hor);

%% Figure 5 - Non Linear IRFs
zz=zeros(1,hor);
figure(5);

for j=1:L
    i=1;
    subplot(2,L,1)
    plot(1:1:hor, zz, 'k-')
    hold on
    plot(1:1:hor, nonlin_res.IRFa.(irfTypes{j})(i,:), 'b--',1:1:hor, nonlin_res.IRFb.(irfTypes{j})(i,:), 'r-o','LineWidth', 1.5)
    axis tight
    ylabel('Government spending')
    xlabel('Horizon (Quarters)')

    i=2;
    subplot(2,L,1+L)
    plot(1:1:hor, zz, 'k-')
    hold on
    plot(1:1:hor, nonlin_res.IRFa.(irfTypes{j})(i,:), 'b--',1:1:hor, nonlin_res.IRFb.(irfTypes{j})(i,:), 'r-o','LineWidth', 1.5)
    axis tight
    ylabel('GDP')
    xlabel('Horizon (Quarters)')

    % add p value plot
end
sgtitle('State Dependent Point Estimates')

%% P value table
est = 'LP';

disp([])

