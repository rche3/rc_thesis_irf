function[outputmultipliers, shock, summaryG, summaryY, x, cum_mult_statea, cum_mult_stateb] = multipliers_irfsmultipliers_irfs_final(lrgdp, lrgov, lrtax, lnlags, snlags, lin_shock, state_shock, FX, fz, experiment, GYsc, hor, fig );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LINEAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=lrgdp(lnlags:end);
g=lrgov(lnlags:end);

Ly=lrgdp(lnlags-1:end-1);
L2y=lrgdp(lnlags-2:end-2);
L3y=lrgdp(lnlags-3:end-3);
L4y=lrgdp(lnlags-4:end-4);

Lg=lrgov(lnlags-1:end-1);
L2g=lrgov(lnlags-2:end-2);
L3g=lrgov(lnlags-3:end-3);
L4g=lrgov(lnlags-4:end-4);

Lt=lrtax(lnlags-1:end-1);
L2t=lrtax(lnlags-2:end-2);
L3t=lrtax(lnlags-3:end-3);
L4t=lrtax(lnlags-4:end-4);

constant= ones(length(y),1);
t=1:1:length(y);
tsq=t.^2;
tcu=t.^3;
tqu=t.^4;

% 
if experiment==3 || experiment==5 || experiment==6
    %Experiment 3, 5, 6
    %[B,BINT,shock] = regress(g,[constant, t', tsq', Lg, L2g, L3g, L4g, Ly, L2y, L3y, L4y, Lt, L2t, L3t, L4t]);
    x=[constant, t', tsq', lin_shock, Lg, L2g, L3g, L4g, Ly, L2y, L3y, L4y, Lt, L2t, L3t, L4t];
    rposl=4;
elseif experiment==0 || experiment==1 || experiment==2
    % %Experiment 0, 1, 2
    %[B,BINT,shock] = regress(g,[constant, Lg, L2g, L3g, L4g, Ly, L2y, L3y, L4y, Lt, L2t, L3t, L4t]);
    x=[constant, lin_shock, Lg, L2g, L3g, L4g, Ly, L2y, L3y, L4y, Lt, L2t, L3t, L4t];
    rposl=2;
elseif experiment==4
    % % %Experiment 4
    %[B,BINT,shock] = regress(g,[constant, t', Lg, L2g, L3g, L4g, Ly, L2y, L3y, L4y, Lt, L2t, L3t, L4t]);
    x=[constant, t', lin_shock, Lg, L2g, L3g, L4g, Ly, L2y, L3y, L4y, Lt, L2t, L3t, L4t];
    rposl=3;
end

[r,nn]=size(x);

for i=1:hor
    yy=y(i:end);
    %yy=(ynew(i:end)- Lynew(1:end-i+1))./Lynew(1:end-i+1);
    results=nwest(yy, x(1:end-i+1,:),i);
    b=results.beta;
    bint(:,1)=b-(results.se*1.96);
    bint(:,2)=b+(results.se*1.96);
    se(:,i)=results.se';
    reglin(:,i)=b;
    confidencelin(:,:,i)=bint;
end
confidencelin=reshape(confidencelin,2*nn,hor);
linGDP=[reglin(rposl,:)];

clear  bint se


for i=1:hor
    yy=g(i:end);
    %yy=(gnew(i:end)- Lgnew(1:end-i+1))./Lynew(1:end-i+1);
    results=nwest(yy, x(1:end-i+1,:),i);
    b=results.beta;
    bint(:,1)=b-(results.se*1.96);
    bint(:,2)=b+(results.se*1.96);
    reg(:,i)=b;
    se(:,i)=results.se';
    confidence(:,:,i)=bint;
end
confidence=reshape(confidence,2*nn,hor);
linGOV=[reg(rposl,:)];


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;%%%%%%%%
% %NON-LINEAR
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear b bint se y g Ly L2y L3y L4y Lg L2g L3g L4g Lt L2t L3t L4t constant t tsq tcu tqu

y=lrgdp(snlags:end);
g=lrgov(snlags:end);

Ly=lrgdp(snlags-1:end-1);
L2y=lrgdp(snlags-2:end-2);
L3y=lrgdp(snlags-3:end-3);
if snlags>=5
    L4y=lrgdp(snlags-4:end-4);
else L4y=[];
end

Lg=lrgov(snlags-1:end-1);
L2g=lrgov(snlags-2:end-2);
L3g=lrgov(snlags-3:end-3);
if snlags>=5
    L4g=lrgov(snlags-4:end-4);
else
    L4g=[];
end

Lt=lrtax(snlags-1:end-1);
L2t=lrtax(snlags-2:end-2);
L3t=lrtax(snlags-3:end-3);
if snlags>=5
    L4t=lrtax(snlags-4:end-4);
else L4t=[];
end

constant= ones(length(y),1);
t=1:1:length(y);
tsq=t.^2;
tcu=t.^3;
tqu=t.^4;

FX=FX(snlags:end,:);
fu=fz(snlags:end);
shock=state_shock;

if experiment==6
    %Experiment 6
    x=[constant, t', tsq', fu, fu.*shock, (1-fu).*shock, fu.*Lg, fu.*L2g, fu.*L3g, fu.*L4g, fu.*Ly, fu.*L2y, fu.*L3y, fu.*L4y,  fu.*Lt, fu.*L2t, fu.*L3t,  fu.*L4t...
        (1-fu).*Lg, (1-fu).*L2g, (1-fu).*L3g, (1-fu).*L4g,  (1-fu).*Ly, (1-fu).*L2y, (1-fu).*L3y, (1-fu).*L4y, (1-fu).*Lt, (1-fu).*L2t, (1-fu).*L3t, (1-fu).*L4t];
    rpos=5;
    [r,nn_upd]=size(x);
elseif experiment==5
    % %Experiment 5
    x=[constant, t', tsq', fu, fu.*shock, (1-fu).*shock, fu.*Lg, fu.*L2g, fu.*L3g, fu.*Ly, fu.*L2y, fu.*L3y,   fu.*Lt, fu.*L2t, fu.*L3t...
        (1-fu).*Lg, (1-fu).*L2g, (1-fu).*L3g, (1-fu).*Ly, (1-fu).*L2y, (1-fu).*L3y, (1-fu).*Lt, (1-fu).*L2t, (1-fu).*L3t];
    [r,nn_upd]=size(x);
    rpos=5;
elseif experiment==4
    % %Experiment 4
    x=[constant, t', fu.*shock, (1-fu).*shock, fu.*Lg, fu.*L2g, fu.*L3g, fu.*Ly, fu.*L2y, fu.*L3y,   fu.*Lt, fu.*L2t, fu.*L3t...
        (1-fu).*Lg, (1-fu).*L2g, (1-fu).*L3g, (1-fu).*Ly, (1-fu).*L2y, (1-fu).*L3y, (1-fu).*Lt, (1-fu).*L2t, (1-fu).*L3t];
    [r,nn_upd]=size(x);
    rpos=3;
elseif experiment==3
    % %Experiment 3
    x=[constant, t', tsq', fu.*shock, (1-fu).*shock, fu.*Lg, fu.*L2g, fu.*L3g, fu.*Ly, fu.*L2y, fu.*L3y,   fu.*Lt, fu.*L2t, fu.*L3t...
        (1-fu).*Lg, (1-fu).*L2g, (1-fu).*L3g, (1-fu).*Ly, (1-fu).*L2y, (1-fu).*L3y, (1-fu).*Lt, (1-fu).*L2t, (1-fu).*L3t];
    [r,nn_upd]=size(x);
    rpos=4;
elseif experiment==2
    % %Experiment 2
    x=[constant,  fu, fu.*shock, (1-fu).*shock, fu.*Lg, fu.*L2g, fu.*L3g, fu.*Ly, fu.*L2y, fu.*L3y,   fu.*Lt, fu.*L2t, fu.*L3t...
        (1-fu).*Lg, (1-fu).*L2g, (1-fu).*L3g, (1-fu).*Ly, (1-fu).*L2y, (1-fu).*L3y, (1-fu).*Lt, (1-fu).*L2t, (1-fu).*L3t];
    [r,nn_upd]=size(x);
    rpos=3;
elseif experiment==1
    % %Experiment 1
    x=[constant, fu.*shock, (1-fu).*shock, fu.*Lg, fu.*L2g, fu.*L3g, fu.*Ly, fu.*L2y, fu.*L3y,   fu.*Lt, fu.*L2t, fu.*L3t...
    (1-fu).*Lg, (1-fu).*L2g, (1-fu).*L3g, (1-fu).*Ly, (1-fu).*L2y, (1-fu).*L3y, (1-fu).*Lt, (1-fu).*L2t, (1-fu).*L3t];
    [r,nn_upd]=size(x);
    rpos=2;
elseif experiment==0
    % %%Experiment 0
    x=[constant, fu.*shock, (1-fu).*shock, fu.*Lg, fu.*L2g, fu.*L3g, fu.*Ly, fu.*L2y, fu.*L3y,   fu.*Lt, fu.*L2t, fu.*L3t...
        (1-fu).*Lg, (1-fu).*L2g, (1-fu).*L3g, (1-fu).*Ly, (1-fu).*L2y, (1-fu).*L3y, (1-fu).*Lt, (1-fu).*L2t, (1-fu).*L3t, FX];
% x=[constant, fu.*shock, (1-fu).*shock, fu.*Lg, fu.*L2g, fu.*L3g, fu.*Ly, fu.*L2y, fu.*L3y,   fu.*Lt, fu.*L2t, fu.*L3t...
%         (1-fu).*Lg, (1-fu).*L2g, (1-fu).*L3g, (1-fu).*Ly, (1-fu).*L2y, (1-fu).*L3y, (1-fu).*Lt, (1-fu).*L2t, (1-fu).*L3t, fu];%, repmat(fu,1,4).*FX, repmat((1-fu),1,4).*FX];
    [r,nn_upd]=size(x);
    rpos=2;
end


for i=1:hor
    yy=y(i:end);
    results=nwest(yy, x(1:end-i+1,:),i);
    b=results.beta;
    bint(:,1)=b-(results.se*1.96);
    bint(:,2)=b+(results.se*1.96);
    se(:,i)=results.se';
    regy(:,i)=b;
    confidencey(:,:,i)=bint;
end
confidencey=reshape(confidencey,2*nn_upd,hor);
stateaGDP=[regy(rpos,:)];
statebGDP=[regy(rpos+1,:)];


clear  b bint se

for i=1:hor
    yy=g(i:end);
    results=nwest(yy, x(1:end-i+1,:),i);
    b=results.beta;
    bint(:,1)=b-(results.se*1.96);
    bint(:,2)=b+(results.se*1.96);
    seg(:,i)=results.se';
    regg(:,i)=b;
    confidenceg(:,:,i)=bint;
end
confidenceg=reshape(confidenceg,2*nn_upd,hor);
stateaGOV=[regg(rpos,:)];
statebGOV=[regg(rpos+1,:)];


outputmultipliers=[linGDP(1,:)./linGOV(1,:); stateaGDP(1,:)./stateaGOV(1,:); statebGDP(1,:)./statebGOV(1,:)];
disp ('Output multipliers for the US')
disp ( 'The multipliers by column are: 2 year integral multiplier, 4 year integral multiplier, 5 year integral multiplier and peak multiplier')
disp('LINEAR ')
[sum(linGDP(1,1:8))./sum(linGOV(1,1:8)), sum(linGDP(1,1:16))./sum(linGOV(1,1:16)), sum(linGDP(1,1:20))./sum(linGOV(1,1:20)),  max(linGDP(1,1:20))./max(linGOV(1,1:20))]* GYsc
disp('HIGH UNEMPLOYMENT/RECESSION state')
[sum(stateaGDP(1,1:8))./sum(stateaGOV(1,1:8)), sum(stateaGDP(1,1:16))./sum(stateaGOV(1,1:16)), sum(stateaGDP(1,1:20))./sum(stateaGOV(1,1:20)), max(stateaGDP(1,1:20))./max(stateaGOV(1,1:20))]* GYsc
disp('LOW UNEMPLOYMENT/NORMAL state')
[sum(statebGDP(1,1:8))./sum(statebGOV(1,1:8)), sum(statebGDP(1,1:16))./sum(statebGOV(1,1:16)),  sum(statebGDP(1,1:20))./sum(statebGOV(1,1:20)), max(statebGDP(1,1:20))./max(statebGOV(1,1:20))]* GYsc

summaryG=[linGOV; stateaGOV; statebGOV]; 
summaryY=[linGDP; stateaGDP; statebGDP]; 



for i=1:hor
    cum_mult_lin(:,i)=  sum(linGDP(1,1:i))./sum(linGOV(1,1:i))*GYsc;
    cum_mult_statea(:,i)= sum(stateaGDP(1,1:i))./sum(stateaGOV(1,1:i))*GYsc;
    cum_mult_stateb(:,i)= sum(statebGDP(1,1:i))./sum(statebGOV(1,1:i))*GYsc;
end
    

figure(11)
plot(1:1:hor, cum_mult_lin, 1:1:hor, cum_mult_statea, '--', 1:1:hor, cum_mult_stateb, '-o', 'LineWidth', 1.5);
legend('Linear', 'Recession', 'Expansion')
xlabel('quarter')
title('Cumulative spending multiplier')
axis tight

confidencelin=confidencelin*GYsc;
reglin=reglin*GYsc;
confidencey=confidencey*GYsc;
regy=regy*GYsc;


zz=zeros(1,hor);
figure(fig)
subplot(2,2,1)
plot(1:1:hor, reg(rposl,:))
hold on
plot(1:1:hor, zz, 'k-')
hold on
grpyat=[(1:1:hor)', confidence(rposl,:)'; (hor:-1:1)' confidence(rposl+nn,hor:-1:1)'];
patch(grpyat(:,1), grpyat(:,2), [0.7 0.7 0.7],'edgecolor', [0.7 0.7 0.7]);%[0.65 0.65 0.65]);
plot(1:1:hor, zz, 'k-')
plot(1:1:hor, reg(rposl,:), 'k','LineWidth', 1.5)
title('Linear: Government spending');
axis tight
subplot(2,2,2)
plot(1:1:hor, reglin(rposl,:))
hold on
plot(1:1:hor, zz, 'k-')
hold on
grpyat=[(1:1:hor)', confidencelin(rposl,:)'; (hor:-1:1)' confidencelin(rposl+nn,hor:-1:1)'];
patch(grpyat(:,1), grpyat(:,2), [0.7 0.7 0.7],'edgecolor', [0.7 0.7 0.7]);%[0.65 0.65 0.65]);
plot(1:1:hor, zz, 'k-')
plot(1:1:hor, reglin(rposl,:), 'k','LineWidth', 1.5)
title('Linear: GDP');
axis tight
subplot(2,2,3)
plot(1:1:hor, regg(rpos,:),'b--')
hold on
grpyat=[(1:1:hor)', confidenceg(rpos,:)'; (hor:-1:1)' confidenceg(rpos+nn_upd,hor:-1:1)'];
patch(grpyat(:,1), grpyat(:,2), [0.7 0.7 0.7],'edgecolor', [0.7 0.7 0.7]);%[0.65 0.65 0.65]);
plot(1:1:hor, regg(rpos,:), 'b--','LineWidth', 1.5)
hold on
plot(1:1:hor, regg(rpos+1,:), 'r-o', 1:1:hor, confidenceg(rpos+1,:), 'r--', 1:1:hor, confidenceg(rpos+1+nn_upd,:), 'r--', 'LineWidth', 1.5);
axis tight
title('State-dependent: Government Spending')
xlabel('quarter')
subplot(2,2,4)
plot(1:1:hor, regy(rpos,:),'b--')
hold on
grpyat=[(1:1:hor)', confidencey(rpos,:)'; (hor:-1:1)' confidencey(rpos+nn_upd,hor:-1:1)'];
patch(grpyat(:,1), grpyat(:,2), [0.7 0.7 0.7],'edgecolor', [0.7 0.7 0.7]);%[0.65 0.65 0.65]);
plot(1:1:hor, regy(rpos,:), 'b--','LineWidth', 1.5)
hold on
plot(1:1:hor, regy(rpos+1,:), 'r-o', 1:1:hor, confidencey(rpos+1,:), 'r--', 1:1:hor, confidencey(rpos+1+nn_upd,:), 'r--', 'LineWidth', 1.5);
axis tight
title('State-dependent: GDP')
xlabel('quarter')
%suptitle('Non-linear case: high unemployment/recession(black), otherwise (red)')
