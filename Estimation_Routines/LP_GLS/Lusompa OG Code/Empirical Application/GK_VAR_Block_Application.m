clear
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code')
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code/Data')
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code/Utilities')

load('Monetary Data Prepped Ramey') %GK Data from Ramey (2016)

tic

ipgrowth=lip(2:end)-lip(1:end-1); %output growth
ipgrowth=[0;ipgrowth];
inflation=lcpi(2:end)-lcpi(1:end-1); %inflation
inflation=[0;inflation];


Y=[ff4_tc gs1 ipgrowth inflation ebp]';

begin=73;

Y=Y(:,begin:end); %data matrix from 1990M1-2012M6


[q, T]=size(Y); 
names=char('mp instrument', '1 year terasury','IP', 'inflation','ebp'); 
Ynames=cellstr(names);
saveY=Y; 

maxh=48; %number of horizons to be estimated not including impact/contemporaneous structural irf (will be maxh+1)
arp=12; %is the lag length
blocksize=10; %the blocksize of the bootstrap
badj=0; %0 for no bias adjustment
nstraps=10001; %number of bootstrap replications

[Structural_IRF,fstat] = var_blockboot_score_iv(Y,arp,nstraps,maxh,blocksize,badj); %SVAR recursive IV Bootstrpap code for structural irfs





%plot for output

i=3;
irf_output=squeeze(Structural_IRF(i,:,:))';
clf;
 plot(1:maxh+1,prctile(irf_output,5),'b:',1:maxh+1,prctile(irf_output,50),'k',1:maxh+1,prctile(irf_output,95),'b:','LineWidth',3)
 xlim([1 maxh+2])
        line([1 maxh+2],[0 0],'color','k','linestyle','--','LineWidth',3); 
        hold off
         box off
        set(gca, 'FontSize', 35)
        title(names(i,:))
        xlabel('Horizon') 
ylabel('Percent')
                grid on


%plot for inflation

i=4;
irf_inflation=squeeze(Structural_IRF(i,:,:))';

figure;
plot(1:maxh+1,prctile(irf_inflation,5),'b:',1:maxh+1,prctile(irf_inflation,50),'k',1:maxh+1,prctile(irf_inflation,95),'b:','LineWidth',3)
 xlim([1 maxh+2])
        line([1 maxh+2],[0 0],'color','k','linestyle','--','LineWidth',3); 
        hold off
         box off
        set(gca, 'FontSize', 35)
        title(names(i,:))
        xlabel('Horizon') 
ylabel('Percent')
                grid on

toc