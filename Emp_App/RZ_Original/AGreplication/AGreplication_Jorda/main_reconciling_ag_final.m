% This code combines part of the ORZ (AEA P&P) code and AG (AEJ) code
clear all
%clc
close all
warning off
%---------------------------------------------------------------------
%           Sample
%---------------------------------------------------------------------
s_desc.startS  = 12;
s_desc.endS    = 252;
s_desc.sample0 = s_desc.startS:1:s_desc.endS;

%---------------------------------------------------------------------
%           			Import data (AG data)
%---------------------------------------------------------------------
% actual macroeconomic time series: output in data_macro
data_macro_series

% actual series for government revenue and spending: output in data0
data_actual_series

% switching series: output in data_SV
data_switching_series

%---------------------------------------------------------------------
%                    create series for VAR (USING AG data)
%---------------------------------------------------------------------
% 
% Core series for the VAR
AG_Y=[data0(s_desc.sample0,5) data0(s_desc.sample0,6) data_macro(s_desc.sample0,5)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment : Choice of Sample size and identification and threshold
% variable and the data used by AG data
% Running our local projection method for 1948.75-2008.75
% Using Perotti-type G shock and F(z) based on moving average of output as
% the threshold variable and exactly the same data used by AG (2012)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lrgov=AG_Y(:,1);
lrtax=AG_Y(:,2);
lrgdp=AG_Y(:,3);

nlags_linear=5;
nlags_state=4;


lin_shock=lrgov(nlags_linear:end);


AGshocks=xlsread('AG_shocks_summary.xlsx',1);
state_shock=AGshocks(1:end,2);%mean
%state_shock=AGshocks(1:end,5);%median
% Switching variable for the VAR
Z=data_SV(s_desc.sample0-1,5); % growth rate of output
% adjust the switching variable to have suffiently long recession periods
% (about 1/4 of the sample).
Z0=Z-0.8;
gamma=-3;

fz=zeros(size(Z));
fz=exp(gamma*Z0(:,1))./(1+exp(gamma*Z0(:,1)));

FX=[data_SV(s_desc.sample0-1,5) data_SV(s_desc.sample0-2,5) data_SV(s_desc.sample0-3,5) data_SV(s_desc.sample0-4,5)];
GYsc=mean(exp(lrgdp)./exp(lrgov));

hor=20;

fig=1;
disp ('Experiment ') % EXPERIMENT 0 replicates the AG scenario
experiment=0
disp('Running Jorda method for 1948.75-2008.75; Using Perotti-type G shock and F(z) based on moving average of output as the threshold ')
disp('variable and exactly the same data used by AG (2012)')
disp('WITH AG SHOCKS')
[outputmultipliers, shock, summaryG, summaryY, x, cum_mult_statea, cum_mult_stateb] = multipliers_irfs_final(lrgdp, lrgov, lrtax, nlags_linear, nlags_state, lin_shock, state_shock, FX, fz, experiment, GYsc, hor, fig);

