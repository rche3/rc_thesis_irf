clear all
close all
clc

%vdata = xlsread('RZ_Data_2014June16.xlsx',1); % US data
%vdata = xlsread('RZ_Data_2014Nov4.xlsx',1); % US data
vdata = xlsread('RZDAT2016_Apr7.xlsx',2); % US data

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
hpunemp=vdata(:,23);


threshold=6.5;

for t=2:size(pdvmil,1)
pdvmily(t,1)= pdvmil(t,1)./ngdp(t-1,1);
end

%     for t=2:size(pdvmil,1)
%         %pdvmily(t,1)= pdvmil(t,1)./(pgdp(t-1,1).*potgdp(t,1));
%         pdvmily(t,1)= pdvmil(t,1)./(pgdp(t-1,1).*potgdp(t-1,1));
%     end

uslrgdp=log(realgdp./totpop);
usrgov=(ngov./totpop./pgdp);
uslrgov=log(ngov./totpop./pgdp);
uslrgov=(uslrgov+ (1-uslrgov(1,1)));
uslrgdp=(uslrgdp+ (1-uslrgdp(1,1)));

   

data=[unemp hpunemp];

rrr=find(unemp>=hpunemp);  % Creating an indicator function with unemployment rate greater than some fixed threshold
dummy=zeros(size(unemp));
dummy(rrr)=1;


ndummy=zeros(size(unemp));
ndummy(rrr)=1;
figure(13)
area(q, ndummy*max(unemp),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
area(q, ndummy*min(unemp),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
plot(q, unemp, q, hpunemp,'k--', 'LineWidth', 1.3)
axis tight




