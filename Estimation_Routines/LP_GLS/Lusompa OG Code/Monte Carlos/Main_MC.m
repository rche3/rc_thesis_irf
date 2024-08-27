%This script reproduces the Monte Carlo Results

clear all
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code')
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code/Monte Carlos')
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code/bias_code/procedures')
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code/Utilities')
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code/Data')


tstart = datetime; 

maxh = 15 ;
maxhLRPY =150 ; %must be greater than or equal to maxh
M=1000;
nstraps = 5001 ;
samplesize=250;
long_auto=16;




%Monte Carlo Results for VAR(1)
coef=[.97 .9 .5];
for i=1:size(coef,2)

[CoverageVAR1, LengthVAR1,BiasVAR1] = MC_VAR1(maxh,coef(i),M,nstraps,samplesize);

v=matlab.lang.makeValidName(['CoverageVAR1_' num2str(coef(i))],'ReplacementStyle','delete');
eval([v '=CoverageVAR1;']);

v=matlab.lang.makeValidName(['LengthVAR1_' num2str(coef(i))],'ReplacementStyle','delete');
eval([v '=LengthVAR1;']);

v=matlab.lang.makeValidName(['BiasVAR1_' num2str(coef(i))],'ReplacementStyle','delete');
eval([v '=BiasVAR1;']);

end



%Monte Carlo Results for arma(1,1)
coef=[0 .25 .5 .75];

for i=1:size(coef,2)

[Coveragearma, Lengtharma,Biasarma] = MC_arma(maxh,coef(i),M,nstraps,samplesize);

v=matlab.lang.makeValidName(['Coveragearma_' num2str(coef(i))],'ReplacementStyle','delete');
eval([v '=Coveragearma;']);

v=matlab.lang.makeValidName(['Lengtharma_' num2str(coef(i))],'ReplacementStyle','delete');
eval([v '=Lengtharma;']);

v=matlab.lang.makeValidName(['Biasarma_' num2str(coef(i))],'ReplacementStyle','delete');
eval([v '=Biasarma;']);

end


%Monte Carlo Results for MA(35) 

[CoverageMA35, LengthMA35,BiasMA35] = MC_MA35(maxh,M,nstraps,samplesize);


%Monte Carlo Results for VARMA(1,1) from Kilian and Kim (2011)

[CoverageVARMA, LengthVARMA,BiasVARMA] = MC_VARMA(maxh,M,nstraps,samplesize);


%Monte Carlo Results for Fiscal VAR

[CoverageFiscalVAR, LengthFiscalVAR, BiasFiscalVAR] = MC_FiscalVAR(maxh,M,nstraps,long_auto);


%Monte Carlo Results for Tech VAR

[CoverageTechVAR, LengthTechVAR, BiasTechVAR] = MC_TechVAR(maxh,M,nstraps,long_auto);


%Monte Carlo Results for Poskit and Yao VARMA

[CoveragePY, LengthPY, BiasPY] = MC_PY_VARMA(maxh,M,nstraps,samplesize);


%Monte Carlo Results for Poskit and Yao VARMA with long run restrictions

[CoveragePYLR, LengthPYLR, BiasPYLR] = MC_PY_VARMA_LR(maxh,maxhLRPY,M,nstraps,samplesize);


%Monte Carlo Results for Sims State Space news model

[CoverageSims, LengthSims,BiasSims] = MC_Sims_SS(maxh,M,nstraps,samplesize);


%Monte Carlo Results for Komunjer and Ng State Space 3 eq NK model

[CoverageKN, LengthKN,BiasKN] = MC_KN_SS(maxh,M,nstraps,samplesize);


tend = datetime; 

tend-tstart

%save results as MCResults.mat. This mat file will be used in to generate
%figures
% fname = '/Users/j1abl01/Desktop/LP GLS Efficiency Code/Monte Carlos'; %uncomment these two lines to save Monte Carlo results
% save(fullfile(fname, 'MCResults.mat'))
