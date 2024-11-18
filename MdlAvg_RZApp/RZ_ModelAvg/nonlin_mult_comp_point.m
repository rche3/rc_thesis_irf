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
irfTypes = {'LP', 'SVAR', 'SVARLPAvg05','SVARLPAvg10', 'SVARLPAvg15', 'SVARLPAvg30'};

L = length(irfTypes);

%% Compute IRFs
% key settings
method = "normal";
p = nlag;
nstraps = 200;


e=1; est=irfTypes{e}; 
% STATE DEPENDENT LP
z=xorig;
x = xorig;
shock = shock;
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
[nlr.IRFa.(est), nlr.IRFb.(est), ~, ~,~, ~, nlr.multa.(est), nlr.multb.(est)]=statelp(data,shock,z,In,hor, ...
                                                                        transformation, clevel, opt, ...
                                                                        bootstrap, method, p, nstraps); 

e=2; est = irfTypes{e};
% SVAR
I = fu;
y = y_tvar;
p = nlag;
rind = [2:3];
sind = 1;
bootstrap = 1;
[nlr.IRFa.(est), nlr.IRFb.(est),~, ~, ~, ~, ~, nlr.multa.(est), nlr.multb.(est)] = stateSVAR(y,I,p,hor,rind,sind, nstraps, bootstrap, clevel, method);

% SVAR / LP MODEL AVG
e=3; est = irfTypes{e};
lambda = ones(1,hor) * 0.05;
[nlr.IRFa.(est), nlr.IRFb.(est), ~,~,~, ~, nlr.multa.(est), nlr.multb.(est)] = stateSVARLP_avg(hor, nlag, nstraps, lambda, method, ...
    data,shock,z, In,transformation, clevel, opt, ...
    y, I, rind, sind);

e=4; est = irfTypes{e};
lambda = ones(1,hor) * 0.1;
[nlr.IRFa.(est), nlr.IRFb.(est), ~,~,~, ~, nlr.multa.(est), nlr.multb.(est)] = stateSVARLP_avg(hor, nlag, nstraps, lambda, method, ...
    data,shock,z, In,transformation, clevel, opt, ...
    y, I, rind, sind);

e=5; est = irfTypes{e};
lambda = ones(1,hor) * 0.15;
[nlr.IRFa.(est), nlr.IRFb.(est), ~,~,~, ~, nlr.multa.(est), nlr.multb.(est)] = stateSVARLP_avg(hor, nlag, nstraps, lambda, method, ...
    data,shock,z, In,transformation, clevel, opt, ...
    y, I, rind, sind);

e=6; est = irfTypes{e};
lambda = ones(1,hor) * 0.3;
[nlr.IRFa.(est), nlr.IRFb.(est), ~,~,~, ~, nlr.multa.(est), nlr.multb.(est)] = stateSVARLP_avg(hor, nlag, nstraps, lambda, method, ...
    data,shock,z, In,transformation, clevel, opt, ...
    y, I, rind, sind);

%% Figure 5 - Non Linear IRFs
close all
zz = zeros(1,hor);
f5 = figure;
ylabels = {'Government Spending', 'GDP'};
% Create empty arrays to hold plot handles
h1 = [];
h2 = [];
% Create arrays to store subplot handles
sp_handles = zeros(2,L);

for j=1:L
    est = irfTypes{j};
    for i=1:2
        sp_handles(i,j) = subplot(2,L,j+L*(i-1));
        hold on
        % Plot IRF lines first
        h1(end+1) = plot(1:1:hor, nlr.IRFa.(est)(i,:), 'b-', 'LineWidth',1.5);
        h2(end+1) = plot(1:1:hor, nlr.IRFb.(est)(i,:), 'r-', 'LineWidth',1.5);
        % Plot zero line last and make it more visible
        plot(1:1:hor, zz, 'k-', 'LineWidth',1, 'HandleVisibility','off')
        axis tight
        ylabel(ylabels{i})
        xlabel('Horizon (Quarters)')
        title(est)
        hold off
    end
end

% Link axes within each row
linkaxes(sp_handles(1,:), 'y')  % Link y-axes for Government Spending row
linkaxes(sp_handles(2,:), 'y')  % Link y-axes for GDP row

sgtitle('State Dependent IRF Estimates')
% Create legend using only one instance of each line type
lgd = legend([h1(1), h2(1)], {'A (ZLB): Point Estimate', 'B (Normal): Point Estimate'}, ...
    'Orientation', 'horizontal');

% Make text larger
lgd.FontSize = 12;  % Adjust this number to make text larger/smaller

% Position the legend at the bottom and make box smaller
% Format is [left bottom width height]
lgd.Position = [0.42, 0.01, 0.18, 0.03];  % Reduce height from 0.05 to 0.03

% Optional: make legend more compact
lgd.ItemTokenSize = [10,10];  % Reduce size of lines in legend
lgd.Box = 'on';  % Remove box around legend (optional)

set(f5, 'Position', [0, 0, 1000, 500])
saveas(f5,'fig/nonlin_mdlavg_pointest.png')

%% Figure 6 - Non Linear Multipliers
close all
zz=zeros(1,hor);
f6 = figure;

% Create empty arrays to hold plot handles
h1 = [];
h2 = [];

for j=1:L
    est = irfTypes{j};
    subplot(1,L,j)
    plot(1:1:hor, zz, 'k-', 'HandleVisibility','off')
    hold on
    h1(end+1) = plot(1:1:hor, nlr.multa.(est), 'b-', 'LineWidth',1.5);
    h2(end+1) = plot(1:1:hor, nlr.multb.(est), 'r-', 'LineWidth',1.5);
    axis tight
    ylabel('Cumulative Multiplier')
    xlabel('Horizon (Quarters)')
    title(est)
end

sgtitle('State Dependent Multiplier Estimates')
% Create legend using only one instance of each line type
lgd = legend([h1(1), h2(1)], {'A (ZLB): Point Estimate', 'B (Normal): Point Estimate'}, ...
    'Orientation', 'horizontal');

% Make text larger
lgd.FontSize = 12;  % Adjust this number to make text larger/smaller

% Position the legend at the bottom and make box smaller
% Format is [left bottom width height]
lgd.Position = [0.42, 0.01, 0.18, 0.03];  % Reduce height from 0.05 to 0.03

% Optional: make legend more compact
lgd.ItemTokenSize = [10,10];  % Reduce size of lines in legend
lgd.Box = 'on';  % Remove box around legend (optional)

set(f6, 'Position', [0, 0, 1200, 500])
saveas(f6,'fig/nonlin_mdlavg_pointest_multiplier.png')

%%

multa_all = [];
multb_all = [];

for j=1:L
    est = irfTypes{j};
    lr_multa = nlr.multa.(est)(end)
    lr_multb = nlr.multb.(est)(end)
    multa_all = [multa_all; lr_multa];
    multb_all = [multb_all; lr_multb];
end

disp('Period A')
multa_all

disp('Period B')
multb_all
