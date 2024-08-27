===========================================================================
DESCRIPTION
This zip file contains MATLAB code for replicating the Monte Carlos in "Local Projections, Autocorrelation, and Efficiency". 
The software used for the code is Matlab. The empirical applications were run in Matlab 2022b. The Monte Carlos were run in 
Matlab 2022a on a cluster with 16 cores and one NVidia V100 GPU.

The VariableDefinitionsandTransformations.txt file describes the data files and transformations used in the empirical VARs in 
the "Monte Carlo" section of the paper as well as the "Empirical Application" in the "Online Appendix".

---------------------------------------------------------------------------
MAIN SCRIPTS

Main_MC.m

    This is the main MATLAB script for Monte Carlos and is located  in the "Monte Carlo" Folder. Code takes approximately 35 
    hours to run using 16 cores and one NVidia V100 GPU. This script requires the distributed computing toolbox, 
    the econometrics toolbox, the optimization toolbox, and the statistics toolbox.

MCResults.mat

    Contains the Monte Carlo results used to produce the figures in the paper. Located in the Monte Carlos folder.

MC_Results_Graphs.m

   This script takes the saved output produced by Main_MC.m (which were saved as MCResults.mat) and produces the figures 
   found in the Monte Carlo section of the paper and in the "Additional Monte Carlo" section of the "Online Appendix". 
   Script is  located in the "Monte Carlo" folder.

GK_LP_GLS_Application.m

    Estimates and creates figures for LP GLS IV application of Gertler and Karadi (2015). Figures can be found in the 
    "Empirical Application" section of the "Online Appendix". Script is located in the "Empirical Application" folder. Code 
    takes approximately 206 seconds to run. Makes use of the statistics toolbox.

GK_LP_OLS_Application.m

    Estimates and creates figures for LP OLS IV application of Gertler and Karadi (2015). Figures can be found in the 
    "Empirical Application" section of the "Online Appendix". Script is located in the "Empirical Application" folder. 
    Code takes approximately 2 seconds to run. Makes use of the statistics toolbox.

GK_VAR_Block_Application.m

    Estimates and creates figures for invertibility robust  SVAR-IV application of Gertler and Karadi (2015). Figures can be 
    found in the "Empirical Application" section of the "Online Appendix". Script is located in the "Empirical Application" 
    folder. Code takes approximately 33 seconds to run. Makes use of the statistics toolbox.
---------------------------------------------------------------------------
SOME IMPORTANT AUXILIARY FUNCTIONS

lp_gls_boot_iv.m

    MATLAB function that runs bootstrap LP GLS IV used in the "Empirical Application" section of the "Online Appendix". Code 
    is annotated and the "How to Section" of the "Online Appendix" walks through how to use the code. Function can be found 
    in the "Utilities" folder.

lp_gls_boot.m            

    MATLAB function that runs bootstrap LP GLS used in the Monte Carlos. Code is annotated and the  "How to Section" of the 
    "Online Appendix" walks through how to use the code. Function can be found in the "Utilities" folder.


===========================================================================