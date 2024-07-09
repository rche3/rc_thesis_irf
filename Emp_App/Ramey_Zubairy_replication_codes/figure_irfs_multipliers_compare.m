clear all
close all
%clc
warning off;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vdata = xlsread('RZDAT.xlsx',2); % US data
time=vdata(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IMPORTANT INPUT (sample window)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if timesample==1
    t1=find(time==1889.75); 
    t2=find(time==2015.75);
elseif timesample==2
    t1=find(time==1947);
    t2=find(time==2015.75);
elseif timesample==3
    t1=find(time==1919);
    t2=find(time==2015.75);
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %DATA
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vdata=vdata(t1:t2,:); 
q=vdata(:,1);
ngov=vdata(:,2);
ngdp=vdata(:,3);
pgdp=vdata(:,4);
totpop=vdata(:,5);
recession=vdata(:,6);
unemp=vdata(:,7);
pdvmil=vdata(:,8);
realgdp=vdata(:,9);
ntax=vdata(:,10);
def=vdata(:,11);
tb3=vdata(:,12);
zlbdummy=vdata(:,13);
potgdp_alt=vdata(:,14);
potgdp=vdata(:,15);

hpunemp=vdata(:,25);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATA TRANSFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating variables of interest
if transformation==1
    for t=2:size(pdvmil,1)
        %pdvmily(t,1)= pdvmil(t,1)./(pgdp(t-1,1).*potgdp(t,1));
        pdvmily(t,1)= pdvmil(t,1)./(pgdp(t-1,1).*potgdp(t-1,1));
    end
elseif transformation==2
    for t=2:size(pdvmil,1)
        pdvmily(t,1)= pdvmil(t,1)./ngdp(t-1,1);
    end
end


rgdp=(realgdp./totpop);
rgov=(ngov./totpop./pgdp);
lrgdp=log(realgdp./totpop);
lrgov=log(ngov./totpop./pgdp);
lrtax=log(ntax./totpop./pgdp);
rgdp_pot=realgdp./potgdp;
rgov_pot=ngov./potgdp./pgdp;
taxy=ntax./ngdp;
defy=def./ngdp;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %IMPORTANT INPUT (threshold for the unemployment rate)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if statechoice==1
    state=unemp;
    threshold=6.5;
elseif statechoice==2
    state=zlbdummy;
    threshold=1;
elseif statechoice==3
    state=unemp;
    threshold=8;
elseif statechoice==4
    state=unemp;
    threshold=hpunemp;
elseif statechoice==5
    state=0.5; % since we want tb3<=0.5
    threshold=tb3;
elseif statechoice==6
    state=recession;
    threshold=1;
end

rrr=find(state>=threshold); 
rrrn=rrr+1; %Since unemployment above threshold in the previous period is the criterion
if rrrn(end)==length(unemp)+1
    rrrn=rrrn(1:end-1);
else
    rrrn=rrrn;
end
fu=zeros(size(unemp));
fu(rrrn)=1; % This is the dummy variable that defines the two states

size(rrr)/size(unemp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if transformation==1 % GK
    if shockchoice==1 %news shock
        if taxchoice==1
            X1=[pdvmily rgdp_pot rgov_pot taxy];
        else
            X1=[pdvmily rgdp_pot rgov_pot];
        end
    else %BP shock
        if taxchoice==1
            X1=[rgdp_pot rgov_pot taxy];
        else
            X1=[rgdp_pot rgov_pot];
        end
    end
elseif transformation==2 % HBR
    if shockchoice==1
        if taxchoice==1
            X1=[pdvmily lrgdp lrgov taxy];
        else
            X1=[pdvmily lrgdp lrgov];
        end
    else
        if taxchoice==1
            X1=[ lrgdp lrgov taxy];
        else
            X1=[ lrgdp lrgov taxy];
        end
    end
end
[rx,cx]=size(X1); 

lags=0:1:nlag;
XLAG1 = lagmatrix(X1,[lags]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SHOCK IDENTIFICATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if shockchoice==1 %news shock
    shock=pdvmily(nlag+1:end);
elseif shockchoice==2  %BP shock
    if transformation==1
        shock=rgov_pot(nlag+1:end);
    elseif transformation==2
        shock=lrgov(nlag+1:end); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRELIMINARIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

constant= ones(length(shock),1);
t=1:1:length(constant);
tsq=t.^2;
tcu=t.^3;
tqu=t.^4;

% %Creating part of the new variables for HBR which are in levels, not logs.
ynew=rgdp(nlag+1:end);
Lynew=rgdp(nlag+1-1:end-1);
gnew=rgov(nlag+1:end);
Lgnew=rgov(nlag+1-1:end-1);


xorig= XLAG1(nlag+1:end, cx+1:end);%only picks lagged values

if transformation==1
    data=[rgov_pot(nlag+1:end), rgdp_pot(nlag+1:end)]; %Y variable
elseif transformation==2
    data=[rgov(nlag+0:end), rgdp(nlag+0:end)];
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LINEAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if trend==4
    x=[constant, t', tsq', tcu', tqu', shock,  xorig]; %complete RHS
    rpos=6; % position of shock
elseif trend==2
    x=[constant, t', tsq', shock,  xorig];
    rpos=4;
elseif trend==0
    x=[constant, shock,  xorig];
    rpos=2;
end

% Generate the standard LP IRFs
[liny, confidencey]=linlp_rz(data,x,hor,rpos,transformation, clevel, opt);


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
% Figure 5 - the IRFs to GDP and Gov. Spending
zz=zeros(1,hor);
i=1;
figure(5)
subplot(2,3,1)
plot(1:1:hor, zz, 'k-')
hold on
plot(1:1:hor, liny(i,:), 'k', 1:1:hor, stateay(i,:), 'b--',1:1:hor, stateby(i,:), 'r-o','LineWidth', 1.5)
axis tight
ylabel('Government spending')
subplot(2,3,2)
grpyat=[(1:1:hor)', confidencey(1,:,i)'; (hor:-1:1)' confidencey(2,hor:-1:1,i)'];
patch(grpyat(:,1), grpyat(:,2), [0.7 0.7 0.7],'edgecolor', [0.7 0.7 0.7]);
hold on
plot(1:1:hor, zz, 'k-')
hold on 
plot(1:1:hor, liny(i,:), 'k','LineWidth', 1.5)
title('Linear')
axis tight
subplot(2,3,3)
grpyat=[(1:1:hor)', confidenceya(1,:,i)'; (hor:-1:1)' confidenceya(2,hor:-1:1,i)'];
patch(grpyat(:,1), grpyat(:,2), [0.7 0.7 0.7],'edgecolor', [0.7 0.7 0.7]);
hold on
plot(1:1:hor, zz, 'k-')
hold on 
plot(1:1:hor, stateay(i,:), 'b--','LineWidth', 1.5)
hold on
plot(1:1:hor, stateby(i,:), 'r-o', 1:1:hor, confidenceyb(1,:,i), 'r--', 1:1:hor, confidenceyb(2,:,i), 'r--', 'LineWidth', 1.5);
axis tight
title('State-dependent')
i=2;
subplot(2,3,4)
plot(1:1:hor, zz, 'k-')
hold on
plot(1:1:hor, liny(i,:), 'k', 1:1:hor, stateay(i,:), 'b--',1:1:hor, stateby(i,:), 'r-o','LineWidth', 1.5)
axis tight
ylabel('GDP')
xlabel('quarter')
subplot(2,3,5)
grpyat=[(1:1:hor)', confidencey(1,:,i)'; (hor:-1:1)' confidencey(2,hor:-1:1,i)'];
patch(grpyat(:,1), grpyat(:,2), [0.7 0.7 0.7],'edgecolor', [0.7 0.7 0.7]);
hold on
plot(1:1:hor, zz, 'k-')
hold on 
plot(1:1:hor, liny(i,:), 'k','LineWidth', 1.5)
xlabel('quarter')
axis tight
subplot(2,3,6)
grpyat=[(1:1:hor)', confidenceya(1,:,i)'; (hor:-1:1)' confidenceya(2,hor:-1:1,i)'];
patch(grpyat(:,1), grpyat(:,2), [0.7 0.7 0.7],'edgecolor', [0.7 0.7 0.7]);
hold on
plot(1:1:hor, zz, 'k-')
hold on 
plot(1:1:hor, stateay(i,:), 'b--','LineWidth', 1.5)
hold on
plot(1:1:hor, stateby(i,:), 'r-o', 1:1:hor, confidenceyb(1,:,i), 'r--', 1:1:hor, confidenceyb(2,:,i), 'r--', 'LineWidth', 1.5);
xlabel('quarter')
axis tight


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Multipliers:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
disp ('Output multipliers for the US')
disp ( 'The multipliers by row are: 2 year integral multiplier, 4 year integral multiplier')
a1=[sum(liny(2,1:8))./sum(liny(1,1:8)); sum(liny(2,1:16))./sum(liny(1,1:16)) ];
a2=[ sum(stateay(2,1:8))./sum(stateay(1,1:8)); sum(stateay(2,1:16))./sum(stateay(1,1:16)) ];
a3=[sum(stateby(2,1:8))./sum(stateby(1,1:8)); sum(stateby(2,1:16))./sum(stateby(1,1:16))];
if statechoice==2
    disp('LINEAR , ZLB state, NORMAL state') 
else
    disp('LINEAR , HIGH UNEMPLOYMENT state, LOW UNEMPLOYMENT state') 
end
    
[a1, a2, a3]

    
cum_mult_lin=  cumsum(liny(2,:))./cumsum(liny(1,:));
cum_mult_statea= cumsum(stateay(2,:))./cumsum(stateay(1,:));
cum_mult_stateb= cumsum(stateby(2,:))./cumsum(stateby(1,:));


if shockchoice==1 %newsy
    std1=xlsread('Multiplier-Standard-Errors.xlsx',2);
    stdlin=std1(1:20,3);
    if statechoice==1 %slack
        std2=xlsread('Multiplier-Standard-Errors.xlsx',2);
        stdb=std2(1:20,5);
        stda=std2(1:20,7);
    elseif statechoice==2 %ZLB
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
multconfl=[cum_mult_lin+clevel*stdlin'; cum_mult_lin-clevel*stdlin'];


figure(6)
subplot(2,1,1)
grpyat=[(1:1:hor)', multconfl(1,:)'; (hor:-1:1)' multconfl(2,hor:-1:1)'];
patch(grpyat(:,1), grpyat(:,2), [0.7 0.7 0.7],'edgecolor', [0.7 0.7 0.7]);
hold on 
plot(1:1:hor, cum_mult_lin, 'k', 'LineWidth', 1.5);
axis tight
title('Linear: cumulative spending multiplier')
subplot(2,1,2)
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



