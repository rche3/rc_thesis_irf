clear
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code')
addpath('/Users/j1abl01/Desktop/LP GLS Efficiency Code/Monte Carlos')
   
load('MCResults')   
fname = '/Users/j1abl01/Desktop/LP GLS Efficiency Code/Figures LP Paper';


    
    maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,Coveragearma_075(:,1)',1:maxh,Coveragearma_075(:,2)',1:maxh,Coveragearma_075(:,3)',1:maxh,Coveragearma_075(:,4)',1:maxh,Coveragearma_075(:,5)',1:maxh,Coveragearma_075(:,6)','LineWidth',4)
xlim([1 maxh])
ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')
        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: ${m=.75}$','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'ARMA .75 Coverage.eps'),'epsc');




    maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,Lengtharma_075(:,1)',1:maxh,Lengtharma_075(:,2)',1:maxh,Lengtharma_075(:,3)',1:maxh,Lengtharma_075(:,4)',1:maxh,Lengtharma_075(:,5)',1:maxh,Lengtharma_075(:,6)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthEast')
title('Average Length: ${m=.75}$','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'ARMA .75 Length.eps'),'epsc');





    maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,Coveragearma_05(:,1)',1:maxh,Coveragearma_05(:,2)',1:maxh,Coveragearma_05(:,3)',1:maxh,Coveragearma_05(:,4)',1:maxh,Coveragearma_05(:,5)',1:maxh,Coveragearma_05(:,6)','LineWidth',4)
xlim([1 maxh])
ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')
        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: ${m=.5}$','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'ARMA .5 Coverage.eps'),'epsc');




    maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,Lengtharma_05(:,1)',1:maxh,Lengtharma_05(:,2)',1:maxh,Lengtharma_05(:,3)',1:maxh,Lengtharma_05(:,4)',1:maxh,Lengtharma_05(:,5)',1:maxh,Lengtharma_05(:,6)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthEast')
title('Average Length: ${m=.5}$','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'ARMA .5 Length.eps'),'epsc');





    maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,Coveragearma_025(:,1)',1:maxh,Coveragearma_025(:,2)',1:maxh,Coveragearma_025(:,3)',1:maxh,Coveragearma_025(:,4)',1:maxh,Coveragearma_025(:,5)',1:maxh,Coveragearma_025(:,6)','LineWidth',4)
xlim([1 maxh])
ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')
        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: ${m=.25}$','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'ARMA .25 Coverage.eps'),'epsc');



    maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,Lengtharma_025(:,1)',1:maxh,Lengtharma_025(:,2)',1:maxh,Lengtharma_025(:,3)',1:maxh,Lengtharma_025(:,4)',1:maxh,Lengtharma_025(:,5)',1:maxh,Lengtharma_025(:,6)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthEast')
title('Average Length: ${m=.25}$','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'ARMA .25 Length.eps'),'epsc');







    maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoverageMA35(:,1)',1:maxh,CoverageMA35(:,2)',1:maxh,CoverageMA35(:,3)',1:maxh,CoverageMA35(:,4)',1:maxh,CoverageMA35(:,5)',1:maxh,CoverageMA35(:,6)','LineWidth',4)
xlim([1 maxh])
% ylim([0 1])
xlabel('Horizon') 
ylabel('Coverage Rate')
        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot','VAR','Location','SouthWest')
title('Coverage: MA(35)','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'MA(35) Coverage.eps'),'epsc');



    maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,LengthMA35(:,1)',1:maxh,LengthMA35(:,2)',1:maxh,LengthMA35(:,3)',1:maxh,LengthMA35(:,4)',1:maxh,LengthMA35(:,5)',1:maxh,LengthMA35(:,6)','LineWidth',4)
xlim([1 maxh])
% ylim([0 1])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot','VAR','Location','SouthWest')
title('Average Length: MA(35)','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'MA(35) Length.eps'),'epsc');








qq=4; %q number of irf per horizon (q*q)   
r=1; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoverageVAR1_097(:,r)',1:maxh,CoverageVAR1_097(:,qq+r)',1:maxh,CoverageVAR1_097(:,(2*qq)+r)',1:maxh,CoverageVAR1_097(:,(3*qq)+r)',1:maxh,CoverageVAR1_097(:,(4*qq)+r)',1:maxh,CoverageVAR1_097(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')
        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: ${A_{11}=.97}$','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'VAR .97 Coverage.eps'),'epsc');



qq=4; %q number of irf per horizon (q*q)   
r=1; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,LengthVAR1_097(:,r)',1:maxh,LengthVAR1_097(:,qq+r)',1:maxh,LengthVAR1_097(:,(2*qq)+r)',1:maxh,LengthVAR1_097(:,(3*qq)+r)',1:maxh,LengthVAR1_097(:,(4*qq)+r)',1:maxh,LengthVAR1_097(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthEast')
title('Average Length: ${A_{11}=.97}$','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'VAR .97 Length.eps'),'eps');




qq=4; %q number of irf per horizon (q*q)   
r=1; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoverageVAR1_09(:,r)',1:maxh,CoverageVAR1_09(:,qq+r)',1:maxh,CoverageVAR1_09(:,(2*qq)+r)',1:maxh,CoverageVAR1_09(:,(3*qq)+r)',1:maxh,CoverageVAR1_09(:,(4*qq)+r)',1:maxh,CoverageVAR1_09(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')
        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: ${A_{11}=.9}$','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'VAR .9 Coverage.eps'),'epsc');



qq=4; %q number of irf per horizon (q*q)   
r=1; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,LengthVAR1_09(:,r)',1:maxh,LengthVAR1_09(:,qq+r)',1:maxh,LengthVAR1_09(:,(2*qq)+r)',1:maxh,LengthVAR1_09(:,(3*qq)+r)',1:maxh,LengthVAR1_09(:,(4*qq)+r)',1:maxh,LengthVAR1_09(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthEast')
title('Average Length: ${A_{11}=.9}$','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'VAR .9 Length.eps'),'epsc');



qq=4; %q number of irf per horizon (q*q)   
r=1; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoverageVAR1_05(:,r)',1:maxh,CoverageVAR1_05(:,qq+r)',1:maxh,CoverageVAR1_05(:,(2*qq)+r)',1:maxh,CoverageVAR1_05(:,(3*qq)+r)',1:maxh,CoverageVAR1_05(:,(4*qq)+r)',1:maxh,CoverageVAR1_05(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')
        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: ${A_{11}=.5}$','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'VAR .5 Coverage.eps'),'epsc');



qq=4; %q number of irf per horizon (q*q)   
r=1; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,LengthVAR1_05(:,r)',1:maxh,LengthVAR1_05(:,qq+r)',1:maxh,LengthVAR1_05(:,(2*qq)+r)',1:maxh,LengthVAR1_05(:,(3*qq)+r)',1:maxh,LengthVAR1_05(:,(4*qq)+r)',1:maxh,LengthVAR1_05(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Average Length: ${A_{11}=.5}$','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'VAR .5 Length.eps'),'epsc');



qq=9; %q number of irf per horizon (q*q)   
r=2; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoverageVARMA(:,r)',1:maxh,CoverageVARMA(:,qq+r)',1:maxh,CoverageVARMA(:,(2*qq)+r)',1:maxh,CoverageVARMA(:,(3*qq)+r)',1:maxh,CoverageVARMA(:,(4*qq)+r)',1:maxh,CoverageVARMA(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')

        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: Inflation Response to Investment','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'VARMA Coverage.eps'),'epsc');




qq=9; %q number of irf per horizon (q*q)   
r=2; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 

    plot(1:maxh,LengthVARMA(:,r)',1:maxh,LengthVARMA(:,qq+r)',1:maxh,LengthVARMA(:,(2*qq)+r)',1:maxh,LengthVARMA(:,(3*qq)+r)',1:maxh,LengthVARMA(:,(4*qq)+r)',1:maxh,LengthVARMA(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','NorthWest')
title('Average Length: Inflation Response to Investment','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'VARMA Length.eps'),'epsc');







qq=4; %q number of irf per horizon (q*q)   
r=1; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoverageFiscalVAR(:,r)',1:maxh,CoverageFiscalVAR(:,qq+r)',1:maxh,CoverageFiscalVAR(:,(2*qq)+r)',1:maxh,CoverageFiscalVAR(:,(3*qq)+r)',1:maxh,CoverageFiscalVAR(:,(4*qq)+r)',1:maxh,CoverageFiscalVAR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
% ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')

        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: Output Response to Output','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'Fiscal Coverage 1.eps'),'epsc');




qq=4; %q number of irf per horizon (q*q)   
r=1; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 

    plot(1:maxh,LengthFiscalVAR(:,r)',1:maxh,LengthFiscalVAR(:,qq+r)',1:maxh,LengthFiscalVAR(:,(2*qq)+r)',1:maxh,LengthFiscalVAR(:,(3*qq)+r)',1:maxh,LengthFiscalVAR(:,(4*qq)+r)',1:maxh,LengthFiscalVAR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Average Length: Output Response to Output','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'Fiscal Length 1.eps'),'epsc');



qq=4; %q number of irf per horizon (q*q)   
r=2; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoverageFiscalVAR(:,r)',1:maxh,CoverageFiscalVAR(:,qq+r)',1:maxh,CoverageFiscalVAR(:,(2*qq)+r)',1:maxh,CoverageFiscalVAR(:,(3*qq)+r)',1:maxh,CoverageFiscalVAR(:,(4*qq)+r)',1:maxh,CoverageFiscalVAR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
% ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')

        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: Spending Response to Output','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'Fiscal Coverage 2.eps'),'epsc');




qq=4; %q number of irf per horizon (q*q)   
r=2; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 

    plot(1:maxh,LengthFiscalVAR(:,r)',1:maxh,LengthFiscalVAR(:,qq+r)',1:maxh,LengthFiscalVAR(:,(2*qq)+r)',1:maxh,LengthFiscalVAR(:,(3*qq)+r)',1:maxh,LengthFiscalVAR(:,(4*qq)+r)',1:maxh,LengthFiscalVAR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Average Length: Spending Response to Output','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'Fiscal Length 2.eps'),'epsc');



qq=4; %q number of irf per horizon (q*q)   
r=3; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoverageFiscalVAR(:,r)',1:maxh,CoverageFiscalVAR(:,qq+r)',1:maxh,CoverageFiscalVAR(:,(2*qq)+r)',1:maxh,CoverageFiscalVAR(:,(3*qq)+r)',1:maxh,CoverageFiscalVAR(:,(4*qq)+r)',1:maxh,CoverageFiscalVAR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
% ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')

        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: Output Response to Spending','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'Fiscal Coverage 3.eps'),'epsc');




qq=4; %q number of irf per horizon (q*q)   
r=3; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 

    plot(1:maxh,LengthFiscalVAR(:,r)',1:maxh,LengthFiscalVAR(:,qq+r)',1:maxh,LengthFiscalVAR(:,(2*qq)+r)',1:maxh,LengthFiscalVAR(:,(3*qq)+r)',1:maxh,LengthFiscalVAR(:,(4*qq)+r)',1:maxh,LengthFiscalVAR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Average Length: Output Response to Spending','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'Fiscal Length 3.eps'),'epsc');






qq=4; %q number of irf per horizon (q*q)   
r=4; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoverageFiscalVAR(:,r)',1:maxh,CoverageFiscalVAR(:,qq+r)',1:maxh,CoverageFiscalVAR(:,(2*qq)+r)',1:maxh,CoverageFiscalVAR(:,(3*qq)+r)',1:maxh,CoverageFiscalVAR(:,(4*qq)+r)',1:maxh,CoverageFiscalVAR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
% ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')

        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: Spending Response to Spending','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'Fiscal Coverage 4.eps'),'epsc');




qq=4; %q number of irf per horizon (q*q)   
r=4; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 

    plot(1:maxh,LengthFiscalVAR(:,r)',1:maxh,LengthFiscalVAR(:,qq+r)',1:maxh,LengthFiscalVAR(:,(2*qq)+r)',1:maxh,LengthFiscalVAR(:,(3*qq)+r)',1:maxh,LengthFiscalVAR(:,(4*qq)+r)',1:maxh,LengthFiscalVAR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Average Length: Spending Response to Spending','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'Fiscal Length 4.eps'),'epsc');









qq=16; %q number of irf per horizon (q*q)   
r=16; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoverageTechVAR(:,r)',1:maxh,CoverageTechVAR(:,qq+r)',1:maxh,CoverageTechVAR(:,(2*qq)+r)',1:maxh,CoverageTechVAR(:,(3*qq)+r)',1:maxh,CoverageTechVAR(:,(4*qq)+r)',1:maxh,CoverageTechVAR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
% ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')

        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: TFP Response to TFP','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'Tech Coverage 1.eps'),'epsc');




qq=16; %q number of irf per horizon (q*q)   
r=16; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 

    plot(1:maxh,LengthTechVAR(:,r)',1:maxh,LengthTechVAR(:,qq+r)',1:maxh,LengthTechVAR(:,(2*qq)+r)',1:maxh,LengthTechVAR(:,(3*qq)+r)',1:maxh,LengthTechVAR(:,(4*qq)+r)',1:maxh,LengthTechVAR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Average Length: TFP Response to TFP','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'Tech Length 1.eps'),'epsc');



qq=16; %q number of irf per horizon (q*q)   
r=1; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoverageTechVAR(:,r)',1:maxh,CoverageTechVAR(:,qq+r)',1:maxh,CoverageTechVAR(:,(2*qq)+r)',1:maxh,CoverageTechVAR(:,(3*qq)+r)',1:maxh,CoverageTechVAR(:,(4*qq)+r)',1:maxh,CoverageTechVAR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
% ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')

        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: Output Response to Output','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'Tech Coverage 2.eps'),'epsc');




qq=16; %q number of irf per horizon (q*q)   
r=1; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 

    plot(1:maxh,LengthTechVAR(:,r)',1:maxh,LengthTechVAR(:,qq+r)',1:maxh,LengthTechVAR(:,(2*qq)+r)',1:maxh,LengthTechVAR(:,(3*qq)+r)',1:maxh,LengthTechVAR(:,(4*qq)+r)',1:maxh,LengthTechVAR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Average Length: Output Response to Output','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'Tech Length 2.eps'),'epsc');



qq=16; %q number of irf per horizon (q*q)   
r=13; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoverageTechVAR(:,r)',1:maxh,CoverageTechVAR(:,qq+r)',1:maxh,CoverageTechVAR(:,(2*qq)+r)',1:maxh,CoverageTechVAR(:,(3*qq)+r)',1:maxh,CoverageTechVAR(:,(4*qq)+r)',1:maxh,CoverageTechVAR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
% ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')

        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: Output Response to TFP','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'Tech Coverage 3.eps'),'epsc');




qq=16; %q number of irf per horizon (q*q)   
r=13; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 

    plot(1:maxh,LengthTechVAR(:,r)',1:maxh,LengthTechVAR(:,qq+r)',1:maxh,LengthTechVAR(:,(2*qq)+r)',1:maxh,LengthTechVAR(:,(3*qq)+r)',1:maxh,LengthTechVAR(:,(4*qq)+r)',1:maxh,LengthTechVAR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Average Length: Output Response to TFP','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'Tech Length 3.eps'),'epsc');






qq=16; %q number of irf per horizon (q*q)   
r=12; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoverageTechVAR(:,r)',1:maxh,CoverageTechVAR(:,qq+r)',1:maxh,CoverageTechVAR(:,(2*qq)+r)',1:maxh,CoverageTechVAR(:,(3*qq)+r)',1:maxh,CoverageTechVAR(:,(4*qq)+r)',1:maxh,CoverageTechVAR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
% ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')

        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: TFP Response to Productivity','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'Tech Coverage 4.eps'),'epsc');




qq=16; %q number of irf per horizon (q*q)   
r=12; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 

    plot(1:maxh,LengthTechVAR(:,r)',1:maxh,LengthTechVAR(:,qq+r)',1:maxh,LengthTechVAR(:,(2*qq)+r)',1:maxh,LengthTechVAR(:,(3*qq)+r)',1:maxh,LengthTechVAR(:,(4*qq)+r)',1:maxh,LengthTechVAR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Average Length: TFP Response to Productivity','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'Tech Length 4.eps'),'epsc');






qq=4; %q number of irf per horizon (q*q)   
r=1; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoveragePY(:,r)',1:maxh,CoveragePY(:,qq+r)',1:maxh,CoveragePY(:,(2*qq)+r)',1:maxh,CoveragePY(:,(3*qq)+r)',1:maxh,CoveragePY(:,(4*qq)+r)',1:maxh,CoveragePY(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
% ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')

        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: Hours Response to Hours','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'PY Coverage.eps'),'epsc');




qq=4; %q number of irf per horizon (q*q)   
r=1; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 

    plot(1:maxh,LengthPY(:,r)',1:maxh,LengthPY(:,qq+r)',1:maxh,LengthPY(:,(2*qq)+r)',1:maxh,LengthPY(:,(3*qq)+r)',1:maxh,LengthPY(:,(4*qq)+r)',1:maxh,LengthPY(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Average Length: Hours Response to Hours','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'PY Length.eps'),'epsc');




qq=4; %q number of irf per horizon (q*q)   
r=3; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoveragePY(:,r)',1:maxh,CoveragePY(:,qq+r)',1:maxh,CoveragePY(:,(2*qq)+r)',1:maxh,CoveragePY(:,(3*qq)+r)',1:maxh,CoveragePY(:,(4*qq)+r)',1:maxh,CoveragePY(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
% ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')

        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: Hours Response to Output','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'PY Coverage 2.eps'),'epsc');




qq=4; %q number of irf per horizon (q*q)   
r=3; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 

    plot(1:maxh,LengthPY(:,r)',1:maxh,LengthPY(:,qq+r)',1:maxh,LengthPY(:,(2*qq)+r)',1:maxh,LengthPY(:,(3*qq)+r)',1:maxh,LengthPY(:,(4*qq)+r)',1:maxh,LengthPY(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Average Length: Hours Response to Output','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'PY Length 2.eps'),'epsc');



qq=4; %q number of irf per horizon (q*q)   
r=1; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoveragePYLR(:,r)',1:maxh,CoveragePYLR(:,qq+r)',1:maxh,CoveragePYLR(:,(2*qq)+r)',1:maxh,CoveragePYLR(:,(3*qq)+r)',1:maxh,CoveragePYLR(:,(4*qq)+r)',1:maxh,CoveragePYLR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
% ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')

        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: Hours Response to Labor Supply Shock','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'PYLR Coverage.eps'),'epsc');




qq=4; %q number of irf per horizon (q*q)   
r=1; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 

    plot(1:maxh,LengthPYLR(:,r)',1:maxh,LengthPYLR(:,qq+r)',1:maxh,LengthPYLR(:,(2*qq)+r)',1:maxh,LengthPYLR(:,(3*qq)+r)',1:maxh,LengthPYLR(:,(4*qq)+r)',1:maxh,LengthPYLR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Average Length: Hours Response to Labor Supply Shock','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'PYLR Length.eps'),'epsc');




qq=4; %q number of irf per horizon (q*q)   
r=3; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoveragePYLR(:,r)',1:maxh,CoveragePYLR(:,qq+r)',1:maxh,CoveragePYLR(:,(2*qq)+r)',1:maxh,CoveragePYLR(:,(3*qq)+r)',1:maxh,CoveragePYLR(:,(4*qq)+r)',1:maxh,CoveragePYLR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
% ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')

        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: Hours Response to Technology Shock','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'PYLR Coverage 2.eps'),'epsc');




qq=4; %q number of irf per horizon (q*q)   
r=3; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 

    plot(1:maxh,LengthPYLR(:,r)',1:maxh,LengthPYLR(:,qq+r)',1:maxh,LengthPYLR(:,(2*qq)+r)',1:maxh,LengthPYLR(:,(3*qq)+r)',1:maxh,LengthPYLR(:,(4*qq)+r)',1:maxh,LengthPYLR(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Average Length: Hours Response to Technology Shock','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'PYLR Length 2.eps'),'epsc');




qq=9; %q number of irf per horizon (q*q)   
r=2; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoverageKN(:,r)',1:maxh,CoverageKN(:,qq+r)',1:maxh,CoverageKN(:,(2*qq)+r)',1:maxh,CoverageKN(:,(3*qq)+r)',1:maxh,CoverageKN(:,(4*qq)+r)',1:maxh,CoverageKN(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
% ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')

        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: Output Response to Interest Rate','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'KN Coverage.eps'),'epsc');




qq=9; %q number of irf per horizon (q*q)   
r=2; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 

    plot(1:maxh,LengthKN(:,r)',1:maxh,LengthKN(:,qq+r)',1:maxh,LengthKN(:,(2*qq)+r)',1:maxh,LengthKN(:,(3*qq)+r)',1:maxh,LengthKN(:,(4*qq)+r)',1:maxh,LengthKN(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Average Length: Output Response to Interest Rate','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'KN Length.eps'),'epsc');



qq=4; %q number of irf per horizon (q*q)   
r=2; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 
    plot(1:maxh,CoverageSims(:,r)',1:maxh,CoverageSims(:,qq+r)',1:maxh,CoverageSims(:,(2*qq)+r)',1:maxh,CoverageSims(:,(3*qq)+r)',1:maxh,CoverageSims(:,(4*qq)+r)',1:maxh,CoverageSims(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
% ylim([0.75 1])
xlabel('Horizon') 
ylabel('Coverage Rate')

        line([1 maxh],[0.95 0.95],'color','k','linestyle','--','LineWidth',3); 

set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Coverage: Output Response to TFP','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'Sims Coverage.eps'),'epsc');




qq=4; %q number of irf per horizon (q*q)   
r=2; %desired response 1-qq
maxh =15; FigH = figure('Position', get(0, 'Screensize')); 

    plot(1:maxh,LengthSims(:,r)',1:maxh,LengthSims(:,qq+r)',1:maxh,LengthSims(:,(2*qq)+r)',1:maxh,LengthSims(:,(3*qq)+r)',1:maxh,LengthSims(:,(4*qq)+r)',1:maxh,LengthSims(:,(5*qq)+r)','LineWidth',4)
xlim([1 maxh])
xlabel('Horizon') 
ylabel('Average Length')
set(gca, 'FontSize', 35)
legend('LP GLS Boot BA','LP GLS Boot','LP GLS','LP OLS','VAR Boot BA','VAR','Location','SouthWest')
title('Average Length: Output Response to TFP','interpreter','latex','FontSize', 35);
saveas(FigH, fullfile(fname, 'Sims Length.eps'),'epsc');
