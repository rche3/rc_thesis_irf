%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vdata = readtable('rzdatnew.csv') % US data
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
% difference the real gdp and gov spending as a multiple of potential gdp
if datafirstdiff == 1
    rgdp_pot = diff(rgdp_pot, 1);
    rgov_pot = diff(rgov_pot, 1);
    pdvmily = pdvmily(2:end);
else
    % pass
end
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

rrr=find(state>=threshold); % returns the indices of the state variable that exceed threshold
rrrn=rrr+1; %Since unemployment above threshold in the previous period is the criterion
if rrrn(end)==length(unemp)+1
    rrrn=rrrn(1:end-1); % truncate if it happens that the final period is a exceeding state
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
XLAG1_lagaug = lagmatrix(X1, [lags, nlag+1]); % lag matrix of "Y" data with an extra lag for lag-augmented LP


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
xorig_lagaug = XLAG1_lagaug(nlag+1:end, cx+1:end);%only picks lagged values

if transformation==1
    data=[rgov_pot(nlag+1:end), rgdp_pot(nlag+1:end)]; % Y variable
elseif transformation==2
    data=[rgov(nlag+0:end), rgdp(nlag+0:end)];
end

%%% create variable for TVAR
if transformation==1
    y_tvar = [pdvmily, rgov_pot, rgdp_pot]; % order as shock, gov, gdp
elseif transformation==2
    y_tvar =[pdvmil, rgov, rgdp];
end