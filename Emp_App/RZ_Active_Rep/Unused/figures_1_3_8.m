clear all
close all
clc

vdata = xlsread('RZDAT.xlsx',2); % US data

time=vdata(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IMPORTANT INPUT (sample window)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1=find(time==1890);
t2=find(time==2015.75);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IMPORTANT INPUT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
tbill=vdata(:,12);
zlbdummy=vdata(:,13);
potgdp_alt=vdata(:,14);
potgdp=vdata(:,15);
nydiscrate=vdata(:,28);

threshold=6.5;

for t=2:size(pdvmil,1)
pdvmily(t,1)= pdvmil(t,1)./(pgdp(t-1,1).*potgdp(t-1,1));
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

%     for t=2:size(pdvmil,1)
%         %pdvmily(t,1)= pdvmil(t,1)./(pgdp(t-1,1).*potgdp(t,1));
%         pdvmily(t,1)= pdvmil(t,1)./(pgdp(t-1,1).*potgdp(t-1,1));
%     end

uslrgdp=log(realgdp./totpop);
usrgov=(ngov./totpop./pgdp);
uslrgov=log(ngov./totpop./pgdp);

nlag=4;
yy=rgov_pot(nlag+1:end);

X1=[rgdp_pot rgov_pot]; 
[rx,cx]=size(X1); 
lags=0:1:nlag;
XLAG1 = lagmatrix(X1,[lags]);
xorig= XLAG1(nlag+1:end, cx+1:end);%only picks lagged values
constant= ones(length(xorig),1);

xx=[constant, xorig];

[B,BINT,bpres] = regress(yy,xx);

bpres=[0;0; 0; 0; bpres];

% 1898q1:  The Spanish-American War starts with the sinking of the USS  Maine.
% 1914q3: WWI starts
% 1939q3:  WWII starts
% 1950q3:  Korean War starts
% 1965q1: Vietnam War starts
% 1980q1: Buildup in response to Soviet invasion of Afghanistan
% 2001q3:  9/11

uslrgov=(uslrgov+ (1-uslrgov(1,1)));
uslrgdp=(uslrgdp+ (1-uslrgdp(1,1)));

figure(1)
subplot(2,1,1)
plot(1890:0.25:2015.75, uslrgov, 'LineWidth', 1)
hold on
vline(1898, 'r--')
vline(1914.5, 'r--')
vline(1939.5, 'r--')
vline(1950.5, 'r--')
vline(1965, 'r--')
vline(1980, 'r--')
vline(2001.5, 'r--')
axis ([1890,2015.75, min(uslrgov), max(uslrgov)+0.25])
title('Log of real per capita government spending ')
subplot(2,1,2)
plot(1890:0.25:2015.75, uslrgdp, 'LineWidth', 1)
hold on
vline(1898, 'r--')
vline(1914.5, 'r--')
vline(1939.5, 'r--')
vline(1950.5, 'r--')
vline(1965, 'r--')
vline(1980, 'r--')
vline(2001.5, 'r--')
axis ([1890,2015.75, min(uslrgdp), max(uslrgdp)+0.25])
title('Log of real per capita GDP ')
    

data=[pdvmily*100 bpres unemp];

rrr=find(unemp>=threshold);  % Creating an indicator function with unemployment rate greater than some fixed threshold
dummy=zeros(size(unemp));
dummy(rrr)=1;


size(rrr)/size(unemp)
%dummy=recession;

varnmz(1,:)='Military news (% of GDP)    ';
varnmz(2,:)='Blanchard-Perotti shock     ';
varnmz(3,:)='Unemployment rate           ';


figure(3)
for i=1:3
    subplot(3,1,i)
    area(q, dummy*max(data(:,i)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
    hold on
    area(q, dummy*min(data(:,i)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
    hold on
    plot(q, data(:,i), 'b', 'LineWidth', 1)
    hold on
    axis ([q(1,1), q(end,1), min(data(:,i)), max(data(:,i))])
    title(varnmz(i,:));
    
end

data=[pdvmily*100 bpres tbill];
varnmz(1,:)='Military news (% of GDP)    ';
varnmz(2,:)='Blanchard-Perotti shock     ';
varnmz(3,:)='3 month Treasury bill rate  ';

%data=[unemp  tbill];
ndummy=[zlbdummy zlbdummy zlbdummy];

figure(8)
for i=1:2
    subplot(3,1,i)
    area(q, ndummy(:,i)*max(data(:,i)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
    hold on
    area(q, ndummy(:,i)*min(data(:,i)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
    hold on
    plot(q, data(:,i), 'b', 'LineWidth', 1)
    hold on
    axis ([q(1,1), q(end,1), min(data(:,i)), max(data(:,i))])
    title(varnmz(i,:));
end
subplot(3,1,3)
i=3;
  area(q, ndummy(:,i)*max(data(:,i)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
    hold on
    area(q, ndummy(:,i)*min(data(:,i)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
    hold on
    plot(q, data(:,i),'b', 'LineWidth', 1)
    hold on
    plot(q, nydiscrate,'b:','LineWidth', 1.5)
    hold on 
    axis ([q(1,1), q(end,1), min(data(:,i)), max(data(:,i))])
    title(varnmz(i,:));


