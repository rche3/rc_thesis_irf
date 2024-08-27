clear;
clc;

%% Introduction
% This is empirical_longhor.m , the main function for replication exercise
% See "explain_longhor.pdf".
% 
% Functions that need to be included to make this code work:
% (1) longhor.m
% (2) longhor1.m
% (3) proc_vb_maq.m
% (4) proc_vb_ma0.m
% (5) proc_vbias.m
% (6) lagmatrix.m
% (7) nwbandwidth.m
% (8) NeweyWest1994.m
% **********************************************************************************************************************************
% 
% Before running,
% in section 2a) give the correct path to open the data set
% in section 2b) input the parameters and choose the suitable series to be xseries yseries (and Wseries in section 3a if needed)
%
% **********************************************************************************************************************************
% Let     Ds(t) = pct change in quarterly exchange rate
%         yy(t) = (Ds(t)+DS(t+1)+...+DS(t+q)/(q+1) = average value of q+1 period change in exchange rate from t-1 to t+q
%         x(t-1) = q+1 period nominal interest differential
% 
% Get bias adjusted estimate of beta in the interest parity regression
% 
%         yy(t) = const. + beta*x(t-1) + eta(t)
% 
% Data are in "empirical_example.xls".  For the UK
%         Ds is called uk_er_3m
%         yy is called  uk_er_6m for q=1, or uk_er_1y for q=4,
%         x is called uk_id_3m, uk_id_6m, etc.
% Names for Japan and Canada are the same, with Canada or jap replacing uk.
% 
% Note the odd dating of yy(t) in the spreadsheet: the entry for, say, 1979:2 and q=1 is [Ds(79:2)+Ds(79:3)]/2
% and thus is realized in 79:3 even though it appears in the 79:2 row in the spreadsheet.  For example,
% for the UK and q=1, uk_er_6m in 1979:2 is 12.08 = (18.81+5.35)/2, where 18.81 and 5.35 are the 79:2 and 79:3
% values of uk_er_3m.
% 
% Needs to be run separately for each currency and horizon in the data set.
%
%**********************************************************************************************************************************
% 
% To apply this file to other data set,
% do a control+F search to search for (@@)
% there are 4 places with (@@), 
% 1) section 2a), change the data set and read in different variables
% 2) section 2b), change parameters accordingly
% 3) section 2b), select the suitable series for VAR estimation and select the relevant xseries yseries (Wseries if needed)
% 4) section 3a), change the position for temp_first and temp_last (first period and last period where all relevant data are available)
% 
% 
% **********************************************************************************************************************************
% 
% As noted above, in the spreadsheet yy(t) is realized in period t+q.  To invoke longhor, the code sets
% 
%         yseries(t+nq) = yy(t)
% 
% so that yseries(t) is realized in period t.
% 
% ********************************************************************************************************************************** 
%
% Structure of the program
% 
% 1. State the purpose of the run
% 2a. Load data
% 2b. Input parameter (equivalent to .inp file in RATS version)
% 3a. declare variables needed by longhor
% 3b. invoke longhor
% 4. compute s.e.
% 5. display results
% require manual input in part 2b.

%% 1. State the purpose of the run
display('Replicate table 6 in the paper using longhor.m');


%% 2a. Load data
% first data point (1) is 1970Q1
% last data point (168) is 2011Q4

% (@@) You can either use Matlab Import Data button or run this code
in_sheet_data = xlsread('example_data.xlsx');

Canada_er_3m=in_sheet_data(:,1);
jap_er_3m=in_sheet_data(:,2);
uk_er_3m=in_sheet_data(:,3);
Canada_id_3m=in_sheet_data(:,4);
jap_id_3m=in_sheet_data(:,5);
uk_id_3m=in_sheet_data(:,6);
Canada_er_6m=in_sheet_data(:,7);
jap_er_6m=in_sheet_data(:,8);
uk_er_6m=in_sheet_data(:,9);
Canada_id_6m=in_sheet_data(:,10);
jap_id_6m=in_sheet_data(:,11);
uk_id_6m=in_sheet_data(:,12);
Canada_er_1y=in_sheet_data(:,13);
jap_er_1y=in_sheet_data(:,14);
uk_er_1y=in_sheet_data(:,15);
Canada_id_1y=in_sheet_data(:,16);
jap_id_1y=in_sheet_data(:,17);
uk_id_1y=in_sheet_data(:,18);
Canada_er_5y=in_sheet_data(:,19);
jap_er_5y=in_sheet_data(:,20);
uk_er_5y=in_sheet_data(:,21);
Canada_id_5y=in_sheet_data(:,22);
jap_id_5y=in_sheet_data(:,23);
uk_id_5y=in_sheet_data(:,24);
Canada_er_10y=in_sheet_data(:,25);
jap_er_10y=in_sheet_data(:,26);
uk_er_10y=in_sheet_data(:,27);
Canada_id_10yc=in_sheet_data(:,28);
jap_id_10yc=in_sheet_data(:,29);
uk_id_10yc=in_sheet_data(:,30);
KKK=in_sheet_data(:,31);
JJJ=in_sheet_data(:,32);
%% 2b. Input parameter (equivalent to .inp file in RATS version)

%(@@)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  parameter input       %%%%%%%%
%%%%%% require input manually %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
first_y=1970;
first_q=1;
start_y=1979;
start_q=2;
end_y=2011; 
end_q=3; % end_q=2 for Canada, end_q=3 for uk or jap
narlag=2;
nvar=2; % this should be nWseries+1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  selection of series   %%%%%%%%
%%%%%% require input manually %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nq=3;  % nq+1=number of quarters, this must be changed when different series are selected
%(@@)
% this is x(t) (interest rate series) in introduction part
xseries=uk_id_1y; % nq+1=number of quarters, for example, if nq=1, then nq+1 = 2quarters = 6m
% choose "uk_id_3m, uk_id_6m, uk_id_1y, uk_id_5y, uk_id_10yc
% jap_id_3m, jap_id_6m, jap_id_1y, jap_id_5y, jap_id_10yc
% Canada_id_3m, Canada_id_6m, Canada_id_1y, Canada_id_5y, Canada_id_10yc"

% this is yy(t) (exchange rate series) in introduction part
yy=uk_er_1y; % nq+1=number of quarters, for example, if nq=1, then nq+1 = 2quarters = 6m
% choose "uk_er_3m, uk_er_6m, uk_er_1y, uk_er_5y, uk_er_10y
% jap_er_3m, jap_er_6m, jap_er_1y, jap_er_5y, jap_er_10y
% Canada_er_3m, Canada_er_6m, Canada_er_1y, Canada_er_5y, Canada_er_10y"

% this is Ds(t) (one period exchange rate) in introduction part
Ds=uk_er_3m;
% choose "uk_er_3m, jap_er_3m, Canada_er_3m"
Wseries=[Ds]; %(@@) additional series can be added in computing moments related to xseries, e.g. Wseries=[Ds Ds2 Ds3];
%if nWseries=0, Wseries will not be used in computing moments
nWseries=1; % (@@) number of columns of Wseries that is used in computing moments related to xseries, 


%% 3a. Construct Input for longhor

nk=1;

%  temp_first and temp_last are the first period and last period where all relevant data are available
%                       example:
%                       t=1         data    N/A     data    N/A    <---(start_y,start_q)
%                       t=2         data    N/A     data    N/A
%                       t=3         data    data    data    data   <---temp_first
%                       t=4         data    data    data    data
%                       t=5         data    data    data    data   <---first is adjusted by temp_first+1+nq+(nk-1) to avoid taking missing data from t=2
%                       t=6         data    data    data    data
%                       t=7         data    data    data    data
%                       t=8         data    data    data    data   <---temp_last
%                       t=9         N/A     data    data    data
%(@@)
temp_first=(start_y-first_y)*4+(start_q-first_q)+1; % 1979Q2, 1970Q1 is first data point in excel
temp_last=(end_y-first_y)*4+(end_q-first_q)+1; %2011Q3
% first and last are adjusted from temp_first and temp_last so no cells with N/A data will be used in estimation
first=temp_first+nq+(nk-1);
last=temp_last;


% As noted above, in the spreadsheet yy(t) is realized in period t+q.  To invoke longhor, the code sets
%          yseries(t+nq) = yy(t)
% so that yseries(t) is realized in period t.

% to create a lagged yseries from yy(t)
y_lag=zeros(size(yy,1),1);
for i=1:size(yy,1)
    if i<=nq
    y_lag(i,1)=NaN;
    else
    y_lag(i,1)=yy(i-nq);
    end
end
yseries=y_lag;

%% 3b invoke longhor
[vbias,betahat,betahat_adj] = longhor(yseries,xseries,Wseries,nWseries,first,last,nq,nk,narlag);

%% 4. compute s.e.

% this is to reproduce regression of interest to get standard error
% this regression is also conducted in longhor/longhor1
con2 = ones(size(yseries,1),1); %constant term
nxlags2 = nq+1:nq+nk; 
nxonly2 = lagmatrix(xseries,nxlags2); %lag of xseries
nX2 = [con2 nxonly2]; % regressor

% these 2 lines are needed to get Tnow (no. of observations for displaying results
[beta_YY,bint1,eta_hat] = regress(yseries(first:last,:),nX2(first:last,:));
Tnow=size(eta_hat,1); 

% compute s.e.
[s0,s1,m] = nwbandwidth(xseries,eta_hat,first,last,nq); % m is the number of lags to include in NW corrected standard errors
[std_NW,t_NW]=NeweyWest1994(yseries(first:last,:),nX2(first:last,:),m); 

%% 5. display results

% display of beta_hat_VAR is hidden in longhor/longhor1, comment out the formatting line if needed 
vbias_over_Tnow=vbias/Tnow;
disp('  vbias     vbias/Tnow');
disp([vbias   vbias_over_Tnow]);
disp('  beta_hat    s.e._NW   beta_hat_adj');
disp([betahat     std_NW(2:nk+1,1)          betahat_adj]);
