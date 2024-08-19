clear all
close all
%clc
%warning off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT SPECIFICATION CHOICES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%timesample=3; % 1 means 1889-2015; 2 means post WWII
timesample=2;
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
    t1=find(time==1914);
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
tbillorig=vdata(:,12);
zlbdummy=vdata(:,13);
potgdp_alt=vdata(:,14);
potgdp=vdata(:,15);

hpunemp=vdata(:,25);
tbillinterp =vdata(:,26);
nydiscrate=vdata(:,27);

if timesample==2
    %tbill=tbillinterp;
    tbill=tbillorig;
   
else
     tbill=tbillorig;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATA TRANSFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rgdp=(realgdp./totpop);
rgov=(ngov./totpop./pgdp);
lrgdp=log(realgdp./totpop);
lrgov=log(ngov./totpop./pgdp);
lrtax=log(ntax./totpop./pgdp);
rgdp_pot=realgdp./potgdp;
rgov_pot=ngov./potgdp./pgdp;
taxy=ntax./ngdp;%ntax./potgdp./pgdp;%
defy=def./ngdp;
lrpot=log(potgdp./totpop);

%CREATING INFLATION
lpgdp=log(pgdp);
inf= (lpgdp(2:end)-lpgdp(1:end-1))*100;

%CREATING YEAR OVER YEAR INFLATION
for i=1:size(pgdp)-4;
    infyoy(i,:)= pgdp(i+4)/pgdp(i)-1;
    %infyoy(i,:)= lpgdp(i+4)-lpgdp(i);
end;
infyoy=infyoy*100;

%CREATING OUTPUT GAP
outgap=lrgdp-lrpot;
outgap=outgap(5:end)*100;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%EXTRA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %TAYLOR RULE RELATED THRESHOLD
% 
% %CREATING A NEW THRESHOLD VARIABLE
% 
tayloruleint=1+ 1.5*infyoy+0.5*outgap;
tbill=tbill(5:end);

figure(11)
subplot(2,1,1)
plot(q(5:end), infyoy, q(5:end), zeros(size(infyoy)));
axis tight
title('inflation (year over year)')
subplot(2,1,2)
plot(q(5:end), outgap, q(5:end), zeros(size(outgap)));
axis tight
title('output gap')

rrr=find(tbill<=tayloruleint);
fu=zeros(size(tbill));
fu(rrr)=1; % This is the dummy variable that defines the two states %1 MEANS STATE OF MONETARY ACCOMODATION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nydiscrate=nydiscrate(5:end);

data=[tayloruleint  tbill infyoy];
dummy=zlbdummy(5:end);
qq=q(5:end);

figure(13)
subplot(3,1,1)
plot(qq, infyoy, qq, zeros(size(infyoy)));
axis tight
title('Inflation (year over year)')
subplot(3,1,2)
plot(qq, outgap, qq, zeros(size(outgap)));
axis tight
title('Output gap')
subplot(3,1,3)
area(qq, dummy*max(data(:,1)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
area(qq, dummy*min(data(:,1)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
plot(qq, tbill, 'b', qq, nydiscrate,'b:', 'LineWidth', 1.5)
hold on
plot(qq, data(:,1), 'r--',  'LineWidth', 1.5)
%legend('T-bill rate', 'Taylor rule int rate')
hold on
plot(qq, zeros(size(outgap)), 'b');
%legend('T-bill rate', 'Taylor rule int rate')
axis ([qq(1,1), qq(end,1), min(data(:,1)), max(data(:,1))])
title('T-bill rate and Taylor rule int rate')



