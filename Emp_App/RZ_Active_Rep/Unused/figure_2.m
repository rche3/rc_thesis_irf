clear all
close all
%clc
warning off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT SPECIFICATION CHOICES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timesample=1; % 1 means 1889-2015; 2 means post WWII

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
    for t=2:size(pdvmil,1)
        %pdvmily(t,1)= pdvmil(t,1)./(pgdp(t-1,1).*potgdp(t,1));
        pdvmily(t,1)= pdvmil(t,1)./(pgdp(t-1,1).*potgdp(t-1,1));
    end



rgdp=(realgdp);
rgov=(ngov./pgdp);
lrgdp=log(realgdp./totpop);
lrgov=log(ngov./totpop./pgdp);
lrtax=log(ntax./totpop./pgdp);
rgdp_pot=realgdp./potgdp;
rgov_pot=ngov./potgdp./pgdp;
taxy=ntax./ngdp *100;
defy=def./ngdp*100;
nypriv=ngdp-ngov;
rypriv=nypriv./pgdp;



x1=nypriv;
x2=ngov;
x3=pdvmily*100;
x4=unemp;
threshold=6.5;
x5=taxy;
x6=defy;
dummy=zlbdummy;

figure(2)
subplot(3,3,1)
i=find(q==1914);
j=find(q==1920);
[ax,p1,p2] = plotyy(q(i:j), x1(i:j), q(i:j), x2(i:j));%,'semilogy','plot');
ylabel(ax(1),'nypriv') % label left y-axis
ylabel(ax(2),'ngov') % label right y-axis
axis(ax(1), [q(i) q(j) min(x1(i:j))-2, max(x1(i:j))+2])
axis(ax(2), [q(i) q(j) min(x2(i:j))-2, max(x2(i:j))+2])
p1.LineStyle = '--';
p1.LineWidth = 2;
p2.LineWidth = 2;
title('Private activity and Government spending')
subplot(3,3,2)
area(q(i:j), dummy(i:j)*max(x3(i:j)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
area(q(i:j), dummy(i:j)*min(x3(i:j)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
plot(q(i:j), x3(i:j), 'LineWidth', 1.5)
title('Military news (% of GDP)')
axis tight
%axis([q(i) q(j) min(x3(i:j))-0.5, max(x3(i:j))+0.5])
subplot(3,3,3)
plot(q(i:j), x4(i:j), 'LineWidth', 1.5)
title('Unemployment rate')
axis tight
subplot(3,3,4)
i=find(q==1938);
j=find(q==1948);
[ax,p1,p2] = plotyy(q(i:j), x1(i:j), q(i:j), x2(i:j));%,'semilogy','plot');
ylabel(ax(1),'nypriv') % label left y-axis
ylabel(ax(2),'ngov') % label right y-axis
axis(ax(1), [q(i) q(j) min(x1(i:j))-2, max(x1(i:j))+2])
axis(ax(2), [q(i) q(j) min(x2(i:j))-2, max(x2(i:j))+2])
p1.LineStyle = '--';
p1.LineWidth = 2;
p2.LineWidth = 2;
subplot(3,3,5)
area(q(i:j), dummy(i:j)*max(x3(i:j)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
area(q(i:j), dummy(i:j)*min(x3(i:j)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
plot(q(i:j), x3(i:j), 'LineWidth', 1.5)
axis tight
%axis([q(i) q(j) min(x3(i:j))-0.5, max(x3(i:j))+0.5])
subplot(3,3,6)
plot(q(i:j), x4(i:j), 'LineWidth', 1.5)
axis tight
subplot(3,3,7)
i=find(q==1948);
j=find(q==1955);
[ax,p1,p2] = plotyy(q(i:j), x1(i:j), q(i:j), x2(i:j));%,'semilogy','plot');
ylabel(ax(1),'nypriv') % label left y-axis
ylabel(ax(2),'ngov') % label right y-axis
axis(ax(1), [q(i) q(j) min(x1(i:j))-2, max(x1(i:j))+2])
axis(ax(2), [q(i) q(j) min(x2(i:j))-2, max(x2(i:j))+2])
p1.LineStyle = '--';
p1.LineWidth = 2;
p2.LineWidth = 2;
subplot(3,3,8)
area(q(i:j), dummy(i:j)*max(x3(i:j)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
area(q(i:j), dummy(i:j)*min(x3(i:j)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
plot(q(i:j), x3(i:j), 'LineWidth', 1.5)
axis tight
%axis([q(i) q(j) min(x3(i:j))-0.5, max(x3(i:j))+0.5])
subplot(3,3,9)
plot(q(i:j), x4(i:j), 'LineWidth', 1.5)
axis tight


figure(100)
subplot(3,4,1)
i=find(q==1914);
j=find(q==1920);
[ax,p1,p2] = plotyy(q(i:j), x1(i:j), q(i:j), x2(i:j));%,'semilogy','plot');
ylabel(ax(1),'nypriv') % label left y-axis
ylabel(ax(2),'ngov') % label right y-axis
axis(ax(1), [q(i) q(j) min(x1(i:j))-2, max(x1(i:j))+2])
axis(ax(2), [q(i) q(j) min(x2(i:j))-2, max(x2(i:j))+2])
p1.LineStyle = '--';
p1.LineWidth = 2;
p2.LineWidth = 2;
%title('Private activity and Government spending')
subplot(3,4,2)
area(q(i:j), dummy(i:j)*max(x3(i:j)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
area(q(i:j), dummy(i:j)*min(x3(i:j)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
plot(q(i:j), x3(i:j), 'LineWidth', 1.5)
title('Military news (% of GDP)')
axis tight
%axis([q(i) q(j) min(x3(i:j))-0.5, max(x3(i:j))+0.5])
subplot(3,4,3)
plot(q(i:j), x4(i:j), 'LineWidth', 1.5)
title('Unemployment rate')
axis tight
subplot(3,4,4)
[ax,p1,p2] = plotyy(q(i:j), x5(i:j), q(i:j), x6(i:j));%,'semilogy','plot');
ylabel(ax(1),'taxy') % label left y-axis
ylabel(ax(2),'defy') % label right y-axis
axis(ax(1), [q(i) q(j) min(x5(i:j))-2, max(x5(i:j))+2])
axis(ax(2), [q(i) q(j) min(x6(i:j))-2, max(x6(i:j))+2])
p1.LineStyle = '--';
p1.LineWidth = 2;
p2.LineWidth = 2;
subplot(3,4,5)
i=find(q==1938);
j=find(q==1948);
[ax,p1,p2] = plotyy(q(i:j), x1(i:j), q(i:j), x2(i:j));%,'semilogy','plot');
ylabel(ax(1),'nypriv') % label left y-axis
ylabel(ax(2),'ngov') % label right y-axis
axis(ax(1), [q(i) q(j) min(x1(i:j))-2, max(x1(i:j))+2])
axis(ax(2), [q(i) q(j) min(x2(i:j))-2, max(x2(i:j))+2])
p1.LineStyle = '--';
p1.LineWidth = 2;
p2.LineWidth = 2;
subplot(3,4,6)
area(q(i:j), dummy(i:j)*max(x3(i:j)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
area(q(i:j), dummy(i:j)*min(x3(i:j)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
plot(q(i:j), x3(i:j), 'LineWidth', 1.5)
axis tight
%axis([q(i) q(j) min(x3(i:j))-0.5, max(x3(i:j))+0.5])
subplot(3,4,7)
plot(q(i:j), x4(i:j), 'LineWidth', 1.5)
axis tight
subplot(3,4,8)
[ax,p1,p2] = plotyy(q(i:j), x5(i:j), q(i:j), x6(i:j));%,'semilogy','plot');
ylabel(ax(1),'taxy') % label left y-axis
ylabel(ax(2),'defy') % label right y-axis
axis(ax(1), [q(i) q(j) min(x5(i:j))-2, max(x5(i:j))+2])
axis(ax(2), [q(i) q(j) min(x6(i:j))-2, max(x6(i:j))+2])
p1.LineStyle = '--';
p1.LineWidth = 2;
p2.LineWidth = 2;
subplot(3,4,9)
i=find(q==1948);
j=find(q==1955);
[ax,p1,p2] = plotyy(q(i:j), x1(i:j), q(i:j), x2(i:j));%,'semilogy','plot');
ylabel(ax(1),'nypriv') % label left y-axis
ylabel(ax(2),'ngov') % label right y-axis
axis(ax(1), [q(i) q(j) min(x1(i:j))-2, max(x1(i:j))+2])
axis(ax(2), [q(i) q(j) min(x2(i:j))-2, max(x2(i:j))+2])
p1.LineStyle = '--';
p1.LineWidth = 2;
p2.LineWidth = 2;
subplot(3,4,10)
area(q(i:j), dummy(i:j)*max(x3(i:j)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
area(q(i:j), dummy(i:j)*min(x3(i:j)),'FaceColor',[0.75 0.75 0.75],'EdgeColor','none')
hold on
plot(q(i:j), x3(i:j), 'LineWidth', 1.5)
axis tight
%axis([q(i) q(j) min(x3(i:j))-0.5, max(x3(i:j))+0.5])
subplot(3,4,11)
plot(q(i:j), x4(i:j), 'LineWidth', 1.5)
axis tight
subplot(3,4,12)
[ax,p1,p2] = plotyy(q(i:j), x5(i:j), q(i:j), x6(i:j));%,'semilogy','plot');
ylabel(ax(1),'taxy') % label left y-axis
ylabel(ax(2),'defy') % label right y-axis
axis(ax(1), [q(i) q(j) min(x5(i:j))-2, max(x5(i:j))+2])
axis(ax(2), [q(i) q(j) min(x6(i:j))-2, max(x6(i:j))+2])
p1.LineStyle = '--';
p1.LineWidth = 2;
p2.LineWidth = 2;


