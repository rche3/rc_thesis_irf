clear all
close all
clc

shockchoice=1; % 1 means news shock; 2 means BP shock

if shockchoice==1
    cifile_rec= xlsread('TVAR_output.xlsx',1);
    cifile_zlb= xlsread('TVAR_output.xlsx',2);
else
    cifile_rec= xlsread('TVAR_output.xlsx',3);
    cifile_zlb= xlsread('TVAR_output.xlsx',4);
end

if shockchoice==1
    recirf=cifile_rec(28:end, 2+3:10);
    nonrecirf=cifile_rec(28:end,13+3:21);
    
    zlbirf=cifile_zlb(28:end-1, 2+3:10);
    nonzlbirf=cifile_zlb(28:end-1,13+3:21);
else
    recirf=cifile_rec(28:end-1, 2:7);
    nonrecirf=cifile_rec(28:end-1,11:16);
    
    zlbirf=cifile_zlb(28:end-1, 2:7);
    nonzlbirf=cifile_zlb(28:end-1,10:15);
end

hor=20;

zz=zeros(20,1);
figure(1)
subplot(2,2,1)
grpyat=[(1:1:hor)', nonrecirf(:,2); (hor:-1:1)' nonrecirf(hor:-1:1,3)];
patch(grpyat(:,1), grpyat(:,2), [0.7 0.7 0.7],'edgecolor', [0.7 0.7 0.7]);
hold on
plot(1:1:hor, zz, 'k-')
hold on 
plot(1:1:hor, nonrecirf(:,1), 'b--','LineWidth', 1.5)
hold on
plot(1:1:hor, recirf(:,1), 'r-o', 1:1:hor, recirf(:,2), 'r--', 1:1:hor,  recirf(:,3), 'r--', 'LineWidth', 1.5);
title('Recession-dependent: Government spending')
axis tight

subplot(2,2,2)
grpyat=[(1:1:hor)', nonrecirf(:,5); (hor:-1:1)' nonrecirf(hor:-1:1,6)];
patch(grpyat(:,1), grpyat(:,2), [0.7 0.7 0.7],'edgecolor', [0.7 0.7 0.7]);
hold on
plot(1:1:hor, zz, 'k-')
hold on 
plot(1:1:hor, nonrecirf(:,4), 'b--','LineWidth', 1.5)
hold on
plot(1:1:hor, recirf(:,4), 'r-o', 1:1:hor, recirf(:,5), 'r--', 1:1:hor,  recirf(:,6), 'r--', 'LineWidth', 1.5);
title('Recession-dependent: GDP')
axis tight


subplot(2,2,3)
grpyat=[(1:1:hor)', nonzlbirf(:,2); (hor:-1:1)' nonzlbirf(hor:-1:1,3)];
patch(grpyat(:,1), grpyat(:,2), [0.7 0.7 0.7],'edgecolor', [0.7 0.7 0.7]);
hold on
plot(1:1:hor, zz, 'k-')
hold on 
plot(1:1:hor, nonzlbirf(:,1), 'b--','LineWidth', 1.5)
hold on
plot(1:1:hor, zlbirf(:,1), 'r-o', 1:1:hor, zlbirf(:,2), 'r--', 1:1:hor,  zlbirf(:,3), 'r--', 'LineWidth', 1.5);
xlabel('quarter')
title('ZLB-dependent: Government spending')
axis tight


subplot(2,2,4)
grpyat=[(1:1:hor)', nonzlbirf(:,5); (hor:-1:1)' nonzlbirf(hor:-1:1,6)];
patch(grpyat(:,1), grpyat(:,2), [0.7 0.7 0.7],'edgecolor', [0.7 0.7 0.7]);
hold on
plot(1:1:hor, zz, 'k-')
hold on 
plot(1:1:hor, nonzlbirf(:,4), 'b--','LineWidth', 1.5)
hold on
plot(1:1:hor, zlbirf(:,4), 'r-o', 1:1:hor, zlbirf(:,5), 'r--', 1:1:hor,  zlbirf(:,6), 'r--', 'LineWidth', 1.5);
xlabel('quarter')
title('ZLB-dependent: GDP')
axis tight