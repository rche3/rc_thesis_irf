clear all
close all
%ZLB figure

ftest=xlsread('Ftests.xlsx',3);%with cumulative G on LHS, difference



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


%Threshold is 23.1085
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ftestcom=xlsread('Ftests.xlsx',5);%MP, modified WWII

hh=ftestcom(:,1);
flincom =ftestcom(:,2);
fexpcom=ftestcom(:,3);
freccom=ftestcom(:,4);
flincomnoww =ftestcom(:,7);
fexpcomnoww=ftestcom(:,8);
freccomnoww=ftestcom(:,9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mm=20;
th=0;
thcom=0;
top= 40 ;
bot=-30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
go=[0, 0.4, 0]; %dark green 
bb=[0.5, 0, 0.9];%purple
or=[1,0.5,0];%brown

figure(10)
subplot(2,3,1)
plot(h(1:mm), flinnews(1:mm), 'LineWidth', 2 , 'Color', bb );
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
subplot(2,3,2)
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
title('ZLB')
axis ([-1,mm, bot,top])
subplot(2,3,3)
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
title('Normal')
axis ([-1,mm, bot,top])

subplot(2,3,4)
plot(h(1:mm), flinnewsnoww(1:mm), 'LineWidth', 2 , 'Color', bb );
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
subplot(2,3,5)
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
subplot(2,3,6)
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






