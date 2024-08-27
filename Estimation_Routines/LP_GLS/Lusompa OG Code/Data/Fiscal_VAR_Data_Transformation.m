clear
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code')
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code/Data')

data=csvread("Fiscal_VAR_Data.csv",1,0 ) 

Date=data(:,1);
ngov=data(:,2);
pgdp=data(:,3);
ngdp=data(:,4);
pop=data(:,5);

logRGDPPC=log(ngdp./(pop.*pgdp));
logRSpendPC=log(ngov./(pop.*pgdp));

RGDPPCGrowth=(logRGDPPC(2:end)-logRGDPPC(1:end-1))*100;
RSpendPCGrowth=(logRSpendPC(2:end)-logRSpendPC(1:end-1))*100;

