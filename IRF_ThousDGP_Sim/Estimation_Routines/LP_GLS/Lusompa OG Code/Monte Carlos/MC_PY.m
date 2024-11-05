clear all
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code')
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code/Monte Carlos')
addpath('/Users/j1abl01/Desktop/Utilities')
addpath('/Users/j1abl01/Desktop/Kolb.VAR.toolbox')
addpath('/Users/j1abl01/Desktop/Boot_Time_Series')
addpath('/Users/j1abl01/Desktop/bias_code/procedures')
addpath('/Users/j1abl01/Desktop/lsbc')

tstart = datetime; 

maxh = 40 ;
M=100;
nstraps = 1001 ;
samplesize=250;
long_auto=16;



% %Monte Carlo Results for Poskit and Yao VARMA
% 
[CoveragePY_cum, LengthPY_cum, BiasPY_cum] = MC_PY_VARMA_Cumulative(maxh,M,nstraps,samplesize);

% %Monte Carlo Results for Poskit and Yao VARMA
% 
% [CoveragePY, LengthPY, BiasPY] = MC_PY_VARMA(maxh,M,nstraps,samplesize);
% 
% rPYVARMA=4;
% 
% %Monte Carlo Results for Poskit and Yao VARMA with long run restrictions
% 
% [CoveragePYLR, LengthPYLR, BiasPYLR] = MC_PY_VARMA_LR(maxh,M,nstraps,samplesize);


tend = datetime; 

tend-tstart
