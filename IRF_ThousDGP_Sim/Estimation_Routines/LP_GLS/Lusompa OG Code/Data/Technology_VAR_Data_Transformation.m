clear
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code')
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code/Data')

data=csvread("Technology_VAR_Data.csv",1,0 ); 
data=data(2:end,:);

Date=data(:,1);
rgdp=data(:,2);
pgdp=data(:,3);
ltfp_util=data(:,4);
pop=data(:,5);
stockp_sh=data(:,6);
tothours=data(:,7);


logrealgdppercapita= log(rgdp./pop);
logrealstockpricespercapita= log(stockp_sh./(pgdp.*pop));
loglaborproductivity =log(rgdp./tothours);

realgdppercapitagrowth =(logrealgdppercapita(2:end)-logrealgdppercapita(1:end-1))*100;
realstockpricespercapitagrowth=(logrealstockpricespercapita(2:end)-logrealstockpricespercapita(1:end-1))*100;
laborproductivitygrowth=(loglaborproductivity(2:end)-loglaborproductivity(1:end-1))*100;
techshockgrowth=(ltfp_util(2:end)-ltfp_util(1:end-1))*100;