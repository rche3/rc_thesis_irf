
%http://www.mathworks.com/matlabcentral/fileexchange/24064-figure-window-generator-fig-m--version-3-1/content/RGB.m
clear all
close all
%Slack figure
ftest=xlsread('Ftests.xlsx',2);%with cumulative G on LHS,  difference

%mpcrit= 23.1085;

h	= ftest(:,1);
flinnews	= ftest(:,2);
fexpnews	= ftest(:,3);
frecnews	= ftest(:,4);
flinbp	= ftest(:,5);
fexpbp	= ftest(:,6);
frecbp= ftest(:,7);
flinnewsnoww	= ftest(:,10);
fexpnewsnoww	= ftest(:,11);
frecnewsnoww	= ftest(:,12);
flinbpnoww	= ftest(:,13);
fexpbpnoww	= ftest(:,14);
frecbpnoww = ftest(:,15);
flinnewspww	= ftest(:,19);
fexpnewspww	= ftest(:,20);
frecnewspww	= ftest(:,21);
flinbppww	= ftest(:,22);
fexpbppww	= ftest(:,23);
frecbppww = ftest(:,24);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ftestcom=xlsread('Ftests_twoinstruments.xlsx',4);%MP
ftestcom=xlsread('Ftests.xlsx',4);%MP, modified WWII, difference

hh=ftestcom(:,1);
flincom =ftestcom(:,2);
fexpcom=ftestcom(:,3);
freccom=ftestcom(:,4);
flincompww =ftestcom(:,7);
fexpcompww=ftestcom(:,8);
freccompww=ftestcom(:,9);
flincomnoww =ftestcom(:,12);
fexpcomnoww=ftestcom(:,13);
freccomnoww=ftestcom(:,14);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mm=20;
th=0;
thcom=0;
top= 40 ;
bot=-30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% colors= {
% 0    0    0    'black';
% 1    0    0    'red';
% 0    1    1    {'cyan','baby blue'};
% 1    0.6  0    'orange';
% 0.5  0.5  1    'light blue';
% 0    1    0    {'green','light green'};
% 0.8  0.5  0    'brown';
% 0.5  0.5  0.5  'dark gray';
% 0.25 0.25 0.9  {'blue','cobalt blue'};
% 1    1    0.6  'cream';
% 0    0.5  0    {'dark green','forest green'};
% 1    0.5  0.5  'peach';
% 1    1    0    'yellow';
% 0    0    0.8  {'dark blue','navy blue'};
% 0.8  0.8  0.8  {'gray','light gray'};
% 0.5  0    0.9  'purple';
% 0.3  0.8  0    'avocado';
% 1    0.5  1    {'magenta','pink'};
% 0    0.8  0.8  {'aqua','turquoise'};
% 0.9  0.75 0    'gold';
% 1    1    1    'white';
% };
go=[0, 0.4, 0]; %dark green 
bb=[0.5, 0, 0.9];%purple
or=[1,0.5,0];%brown

figure(4)
subplot(3,3,1)
plot(h(1:mm), flinnews(1:mm), 'LineWidth',2, 'Color', bb );
hold on
plot(h(1:mm), flinbp(1:mm), '--',   'LineWidth', 2 , 'Color', go);%, 'MarkerEdgeColor','r', 'MarkerFaceColor','r','MarkerSize',5);
hold on
plot(h(1:mm), flincom(1:mm), '-*', 'LineWidth', 1.5, 'Color', or);
hold on
plot(h(1:mm), th*ones(mm,1), 'k--', 'LineWidth', 1)
hold on
plot(h(1:mm), thcom*ones(mm,1), 'k-.', 'LineWidth', 1)
%xlabel('h')
ylabel('Full Sample')
title('Linear')
axis ([-1,mm, bot,top])
subplot(3,3,2)
plot(h(1:mm), frecnews(1:mm), 'LineWidth', 2 , 'Color', bb );
hold on
plot(h(1:mm), frecbp(1:mm), '--',   'LineWidth', 2 , 'Color', go);%, 'MarkerEdgeColor','r', 'MarkerFaceColor','r','MarkerSize',5);
hold on
plot(h(1:mm), freccom(1:mm), '-*', 'LineWidth', 1.5, 'Color', or);
hold on
plot(h(1:mm), th*ones(mm,1), 'k--', 'LineWidth', 1)
hold on
plot(h(1:mm), thcom*ones(mm,1), 'k-.', 'LineWidth', 1)
%xlabel('h')
title('High Unemployment')
axis ([-1,mm, bot,top])
subplot(3,3,3)
plot(h(1:mm), fexpnews(1:mm), 'LineWidth', 2 , 'Color', bb );
hold on
plot(h(1:mm), fexpbp(1:mm), '--',   'LineWidth', 2 , 'Color', go );%, 'MarkerEdgeColor','r', 'MarkerFaceColor','r','MarkerSize',5);
hold on
plot(h(1:mm), fexpcom(1:mm), '-*', 'LineWidth', 1.5, 'Color', or);
hold on
plot(h(1:mm), th*ones(mm,1), 'k--', 'LineWidth', 1)
hold on
plot(h(1:mm), thcom*ones(mm,1), 'k-.', 'LineWidth', 1)
%xlabel('h')
title('Low Unemployment')
axis ([-1,mm, bot,top])

subplot(3,3,4)
plot(h(1:mm), flinnewspww(1:mm), 'LineWidth', 2 , 'Color', bb );
hold on
plot(h(1:mm), flinbppww(1:mm), '--',   'LineWidth', 2 , 'Color', go);%, 'MarkerEdgeColor','r', 'MarkerFaceColor','r','MarkerSize',5);
hold on
plot(h(1:mm), flincompww(1:mm), '-*', 'LineWidth', 1.5, 'Color', or);
hold on
plot(h(1:mm), th*ones(mm,1), 'k--', 'LineWidth', 1)
hold on
plot(h(1:mm), thcom*ones(mm,1), 'k-.', 'LineWidth', 1)
%xlabel('h')
ylabel('Post-WWII')
%title('Linear')
axis ([-1,mm, bot,top])
subplot(3,3,5)
plot(h(1:mm), frecnewspww(1:mm), 'LineWidth', 2 , 'Color', bb );
hold on
plot(h(1:mm), frecbppww(1:mm),'--',   'LineWidth', 2 , 'Color', go );%, 'MarkerEdgeColor','r', 'MarkerFaceColor','r','MarkerSize',5);
hold on
plot(h(1:mm), freccompww(1:mm), '-*', 'LineWidth', 1.5, 'Color', or);
hold on
plot(h(1:mm), th*ones(mm,1), 'k--', 'LineWidth', 1)
hold on
plot(h(1:mm), thcom*ones(mm,1), 'k-.', 'LineWidth', 1)
%xlabel('h')
%title('Recession')
axis ([-1,mm, bot,top])
subplot(3,3,6)
plot(h(1:mm), fexpnewspww(1:mm), 'LineWidth', 2 , 'Color', bb );
hold on
plot(h(1:mm), fexpbppww(1:mm), '--',   'LineWidth', 2 , 'Color', go);%, 'MarkerEdgeColor','r', 'MarkerFaceColor','r','MarkerSize',5);
hold on
plot(h(1:mm), fexpcompww(1:mm), '-*', 'LineWidth', 1.5, 'Color', or);
hold on
plot(h(1:mm), th*ones(mm,1), 'k--', 'LineWidth', 1)
hold on
plot(h(1:mm), thcom*ones(mm,1), 'k-.', 'LineWidth', 1)
%xlabel('h')
%title('Expansion')
axis ([-1,mm, bot,top])

subplot(3,3,7)
plot(h(1:mm), flinnewsnoww(1:mm), 'LineWidth',2, 'Color', bb );
hold on
plot(h(1:mm), flinbpnoww(1:mm), '--',   'LineWidth', 2 , 'Color', go );%, 'MarkerEdgeColor','r', 'MarkerFaceColor','r','MarkerSize',5);
hold on
plot(h(1:mm), flincomnoww(1:mm), '-*', 'LineWidth', 1.5, 'Color', or);
hold on
plot(h(1:mm), th*ones(mm,1), 'k--', 'LineWidth', 1)
hold on
plot(h(1:mm), thcom*ones(mm,1), 'k-.', 'LineWidth', 1)
xlabel('h')
ylabel('Excluding WWII')
%title('Linear')
axis ([-1,mm, bot,top])
subplot(3,3,8)
plot(h(1:mm), frecnewsnoww(1:mm), 'LineWidth', 2 , 'Color', bb);
hold on
plot(h(1:mm), frecbpnoww(1:mm), '--',   'LineWidth', 2 , 'Color', go );%, 'MarkerEdgeColor','r', 'MarkerFaceColor','r','MarkerSize',5);
hold on
plot(h(1:mm), freccomnoww(1:mm), '-*', 'LineWidth', 1.5, 'Color', or);
hold on
plot(h(1:mm), th*ones(mm,1), 'k--', 'LineWidth', 1)
hold on
plot(h(1:mm), thcom*ones(mm,1), 'k-.', 'LineWidth', 1)
xlabel('h')
%title('Recession')
axis ([-1,mm, bot,top])
subplot(3,3,9)
plot(h(1:mm), fexpnewsnoww(1:mm), 'LineWidth', 2 , 'Color', bb);
hold on
plot(h(1:mm), fexpbpnoww(1:mm), '--',   'LineWidth', 2 , 'Color', go );%, 'MarkerEdgeColor','r', 'MarkerFaceColor','r','MarkerSize',5);
hold on
plot(h(1:mm), fexpcomnoww(1:mm), '-*', 'LineWidth', 1.5, 'Color', or);
hold on
plot(h(1:mm), th*ones(mm,1), 'k--', 'LineWidth', 1)
hold on
plot(h(1:mm), thcom*ones(mm,1), 'k-.', 'LineWidth', 1)
xlabel('h')
legend('Military news shock', 'Blanchard-Perotti shock', 'Combined');
%title('Expansion')
axis ([-1,mm, bot,top])

