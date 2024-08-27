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


%calculating structural irfs and conf intervals for output with 90% conf
%intervals
i=3;
irf_output=[];

for h=1:maxh+1

        p=2+arp*q; 
F=zeros(q*arp,T-arp-h+1); 
for j=1:arp 
    F((j-1)*q+(1:q),:)=saveY(:,(arp-j+1):(T-h+1-j)); 
end

saveF=[ones(1,T-arp-h+1);saveY(2,arp+1:end-h+1);F];
yeff=saveY(:,arp+h:end)'; 

[betahat,~,se_beta,~,~,q_df] = har_iv(yeff(:,i),saveF',Y(1,arp+1:end-h+1)',0);

irf_output= [irf_output; betahat(2)-tinv(.95,q_df)*se_beta(2) betahat(2) betahat(2)+tinv(.95,q_df)*se_beta(2)];

end


%plot for output
clf

 plot(1:maxh+1,irf_output(:,1),'b:',1:maxh+1,irf_output(:,2),'k',1:maxh+1,irf_output(:,3),'b:','LineWidth',3)
 xlim([1 maxh+2])
        line([1 maxh+2],[0 0],'color','k','linestyle','--','LineWidth',3); 
        hold off
         box off
        set(gca, 'FontSize', 35)
        title(names(i,:))
        xlabel('Horizon') 
ylabel('Percent')
                grid on


%calculating structural irfs and conf intervals for inflation with 90% conf
%intervals
                i=4;
irf_inflation=[];

for h=1:maxh+1

        p=2+arp*q; 
F=zeros(q*arp,T-arp-h+1); 
for j=1:arp 
    F((j-1)*q+(1:q),:)=saveY(:,(arp-j+1):(T-h+1-j)); 
end

saveF=[ones(1,T-arp-h+1);saveY(2,arp+1:end-h+1);F];
yeff=saveY(:,arp+h:end)'; 

[betahat,~,se_beta,~,~,q_df] = har_iv(yeff(:,i),saveF',Y(1,arp+1:end-h+1)',0);

irf_inflation= [irf_inflation; betahat(2)-tinv(.95,q_df)*se_beta(2) betahat(2) betahat(2)+tinv(.95,q_df)*se_beta(2)];

end


%plot for inflation
figure;
 plot(1:maxh+1,irf_inflation(:,1),'b:',1:maxh+1,irf_inflation(:,2),'k',1:maxh+1,irf_inflation(:,3),'b:','LineWidth',3)
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
