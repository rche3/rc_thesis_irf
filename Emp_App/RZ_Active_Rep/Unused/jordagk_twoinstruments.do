**** JORDAGK.DO

*** May 22, 2016
*** This program estimates IRFs and both 2-step and 1-step multipliers using the Jorda method

*** Uses Gordon-Krenn (GK) transformation (i.e. variables are divided by potential)

***
*** Requires:
***     rzdat.xlsx, updated April 7, 2016
********************************************************************************

 #delimit;

drop _all;
clear all;

set more 1;
set matsize 800;

capture log close;
log using jordagkirfs_results.log, replace;


/*******************************************************************************
  SET PARAMETERS THAT GOVERN SPECIFICATION
*******************************************************************************/

local sample = 1;  /*1 = full sample, 2 = post-WWII */

local omit wwii;  /*either nomit(don't omit subsample) or wwii (omit WWII) */

local state zlb;  /* slack or zlb or recession */

local shock newsy; /* shock identification: either newsy or bp */

local p = 4; /*number of lags of control variables*/

local trends = 0; /*0 = no trends, 1 = trends */

local tax = 0; /*0 = exclude taxes, 1 = include taxes */


*******************************************************************************;
** RAW DATA IMPORTATION AND DATA SETUP;
*******************************************************************************;

import excel rzdat.xlsx, sheet("rzdat") firstrow;

*insheet using rzdatnew.csv;

drop if quarter<1889;

gen qdate = q(1889q1) + _n-1;
tsset qdate, q;

/* World War II rationing sample.  Start in 1941q3 because of constraints on
     auto production.  Most rationing ended in Aug/Sept 1945, a few items in
	 Nov/Dec 1945 */

gen wwii = quarter>=1941.5 & quarter<1946;

gen nomit = 0;  /* indicator for no omit */


*** DEFINE QUARTIC TREND;

gen t = _n;
gen t2 = t^2;
gen t3 = t^3;
gen t4 = t^4;

*** DEFINE STATE VARIABLE;

gen slack = unemp >= 6.5;  /* unemployment state with fixed threshold */
gen slack8 = unemp>=8;
gen slackhp = unemp>=hpunemp_split;  /* unemployment state with hp threshold */
gen zlb5 =tbill<=0.5;
gen zlb = zlb_dummy;  /*  zlb state */


*** NORMALIZATION;

/* choice of potential GDP for normalization:

   rgdp_potCBO (cubic trend early, CBO late) or rgdp_pott6 (6th degree for full),
   both fitted excluding Great Depression, WWII:  quarter>=1930 & quarter<1947*/

local ynorm rgdp_pott6; /* rgdp_pott6 or rgdp_potcbo */

* BASIC VARIABLES;

gen newsy = news/(L.`ynorm'*L.pgdp);
gen rgov = ngov/pgdp;
gen rtax = nfedcurrreceipts_nipa/pgdp;
gen taxy = nfedcurrreceipts_nipa/ngdp;
gen debty = pubfeddebt_treas/L.ngdp;
gen lpgdp = ln(pgdp);

gen infl = 400*D.lpgdp;

* normalize variables and shorten names;

gen y = rgdp/`ynorm';
gen g = rgov/`ynorm';
 
gen bp = g; /* Blanchard-Perotti shock is just orthogonalized current g */

*******************************************************************************;
** CUMULATIVE VARIABLES;
*******************************************************************************;

gen cumuly = 0;
gen cumulg = 0;
 
forvalues i = 0/20 {;

   gen f`i'cumuly = F`i'.y + cumuly;
   gen f`i'cumulg = F`i'.g + cumulg;
   
   gen recf`i'cumulg = f`i'cumulg*L.`state';
   gen expf`i'cumulg = f`i'cumulg*(1-L.`state');
   
   replace cumuly = f`i'cumuly;
   replace cumulg = f`i'cumulg;
   
};


*******************************************************************************;
**  INTERACTION OF SHOCKS WITH STATE;
*******************************************************************************;

 foreach var in newsy bp { ;
 
   gen rec0`var' = `var'*L.`state';
   gen exp0`var' = `var'*(1-L.`state');
 
 };

*******************************************************************************;
** CREATE LISTS;
*******************************************************************************;

   if `sample'==1 {;

         gen h = t - 1;  /* h is the horizon for the irfs */
         global trendlist t t2 t3 t4;
   };

    else {;
         drop if quarter<1947;
         gen h = t - 232 - 1;
         global trendlist t t2;
     };
	 
forvalues i = 1/`p' {; 

  foreach var in newsy y g taxy debty {;

    gen rec`var'`i' = L`i'.`var'*L.`state';
    gen exp`var'`i' = L`i'.`var'*(1-L.`state');
 
  };
};

  if `trends'==0 {;
  
    if `tax'==0 {;
  
      global newsylinxlist L(1/`p').newsy L(1/`p').y L(1/`p').g ;
      global bplinxlist L(1/`p').y L(1/`p').g ;
	  global newsynlxlist L.`state' recnewsy? expnewsy? recy? expy? recg? expg? ;
      global bpnlxlist L.`state' recy? expy? recg? expg? ;
	
	};
	
	else {;
	
	  global newsylinxlist L(1/`p').newsy L(1/`p').y L(1/`p').g L(1/`p').taxy;
      global bplinxlist L(1/`p').y L(1/`p').g L(1/`p').taxy;
	  global newsynlxlist L.`state' recnewsy? expnewsy? recy? expy? recg? expg? rectaxy? exptaxy?;
      global bpnlxlist L.`state' recy? expy? recg? expg? rectaxy? exptaxy?;
	
    };
  };
  
  else {;
  
    if `tax'==0 {;
	
      global newsylinxlist L(1/`p').newsy L(1/`p').y L(1/`p').g $trendlist;
      global bplinxlist L(1/`p').y L(1/`p').g $trendlist;
      global newsynlxlist L.`state' recnewsy? expnewsy? recy? expy? recg? expg? $trendlist;
      global bpnlxlist L.`state' recy? expy? recg? expg? $trendlist;
	
	};
	
	else {;
	
	global newsylinxlist L(1/`p').newsy L(1/`p').y L(1/`p').g L(1/`p').taxy $trendlist;
    global bplinxlist L(1/`p').y L(1/`p').g L(1/`p').taxy $trendlist;
	global newsynlxlist L.`state' recnewsy? expnewsy? recy? expy? recg? expg? rectaxy? exptaxy? $trendlist;
    global bpnlxlist L.`state' recy? expy? recg? expg? rectaxy? exptaxy? $trendlist;
	
    };
	
};


global newsylinshock newsy;
global newsynlshock rec0newsy exp0newsy;

global bplinshock bp;
global bpnlshock rec0bp exp0bp;


** INITIALIZE SUM OF EFFECTS TO 0 AND PARAMETERS SERIES TO MISSING;

gen sumliny = 0; gen sumling = 0;
gen sumexpy = 0; gen sumexpg = 0;
gen sumrecy = 0; gen sumrecg = 0;

foreach var in bylin byexp byrec bglin bgexp bgrec up95bylin up95byexp up95byrec up95bglin up95bgexp up95bgrec
  lo95bylin lo95byexp lo95byrec lo95bglin lo95bgexp lo95bgrec seylin seyexp seyrec seglin segexp segrec
  multlin multexp multrec {;
  
  quietly gen `var' = .;
  
}; 


*****************************************************************************;
drop sey*;


foreach var in multlin1 multexp1 multrec1 Fkplin Fkpexp Fkprec seylin seyexp seyrec ptestdiff multexpb multrecb Fefflin Feffexp Feffrec Feffboth Critlin Critexp Critrec Critboth arp Fdifflin Fdiffexp Fdiffrec{;
  
  quietly gen `var' = .;
  
}; 
*******************************************************************************;
** ESTIMATION OF CUMULATIVE;
*******************************************************************************;

forvalues i = 1/20 {; 

  ivreg2 f`i'cumuly (f`i'cumulg = $`shock'linshock bp) $`shock'linxlist if L`p'.`omit'==0 & F`i'.`omit'==0 & `omit'==0, robust bw(auto);
  gen Fkplinh`i'= e(widstat); /* Kleibergen-Paap rk Wald F statistic*/
  gen multlinh`i' = _b[f`i'cumulg];
  gen seylinh`i' = _se[f`i'cumulg]; /* HAC robust standard error*/
  
  weakivtest;
  gen Fefflinh`i' = r(F_eff);
  gen Critlinh`i'= r(c_TSLS_10);
  gen Fdifflinh`i'= r(F_eff)-r(c_TSLS_10);
  
  ivreg2 f`i'cumuly (expf`i'cumulg = exp0`shock' exp0bp) $`shock'nlxlist if L`p'.`omit'==0 & F`i'.`omit'==0 & `omit'==0, robust bw(auto);
  gen Fkpexph`i'= e(widstat);
  gen multexph`i' = _b[expf`i'cumulg];
  gen seyexph`i' = _se[expf`i'cumulg];
  
  weakivtest;
  gen Feffexph`i' = r(F_eff);
  gen Critexph`i'= r(c_TSLS_10);
  gen Fdiffexph`i'= r(F_eff)-r(c_TSLS_10);
  
  ivreg2 f`i'cumuly (recf`i'cumulg = rec0`shock' rec0bp) $`shock'nlxlist if L`p'.`omit'==0 & F`i'.`omit'==0 & `omit'==0, robust bw(auto);
  gen Fkprech`i'= e(widstat);  
  gen multrech`i' = _b[recf`i'cumulg];
  gen seyrech`i' = _se[recf`i'cumulg];
  
  weakivtest;
  gen Feffrech`i' = r(F_eff);
  gen Critrech`i'= r(c_TSLS_10);
  gen Fdiffrech`i'= r(F_eff)-r(c_TSLS_10);
  
  ivreg2 f`i'cumuly (expf`i'cumulg recf`i'cumulg = exp0`shock' rec0`shock' exp0bp rec0bp) $`shock'nlxlist if L`p'.`omit'==0 & F`i'.`omit'==0 & `omit'==0, robust bw(auto);
  test expf`i'cumulg=recf`i'cumulg;	
  gen ptestdiffh`i' = r(p);
  gen multexpbh`i' = _b[expf`i'cumulg];
  gen multrecbh`i' = _b[recf`i'cumulg];
  
  weakiv ivreg2 f`i'cumuly (f`i'cumulg recf`i'cumulg = $`shock'linshock rec0`shock' exp0bp rec0bp) $`shock'nlxlist if L`p'.`omit'==0 & F`i'.`omit'==0 & `omit'==0, robust bw(auto) strong(f`i'cumulg);
  /*test recf`i'cumulg;*/
  gen arph`i' = e(ar_p);
 
 
    
 foreach var in multlin multexp multrec {;
  
    quietly replace `var'1 = `var'h`i' if h==`i';
	
  };
  
  foreach var in seylin seyexp seyrec ptestdiff multexpb multrecb arp{;
  
    quietly replace `var' = `var'h`i' if h==`i';
	
  };
  
 foreach var in  Fkplin Fkpexp Fkprec Fefflin Feffexp Feffrec Critlin Critexp Critrec{;
  
    quietly replace `var' = `var'h`i' if h==`i';
		
  };
  
  foreach var in  Fdifflin Fdiffexp Fdiffrec {;
  
    quietly replace `var' = `var'h`i' if h==`i';
	quietly replace `var' = 30 if `var'>30;
	
  };
  
};


display as text "First stage F-statistic (Kleibergen-Paap rk Wald F statistic): Linear, Expansion, Recession";


list h Fkplin Fkpexp Fkprec if h<=20;
outsheet h Fdifflin Fdiffexp Fdiffrec using junkfdiff.csv if h<=20, comma replace ;
outsheet h Fefflin Feffexp Feffrec using junkmp.csv if h<=20, comma replace ;

display as text "Montiel-Pflueger robust weak instrument test and TSLS critical values for tau=10%";

list h Fefflin Critlin Feffexp Critexp Feffrec Critrec Feffboth Critboth if h<=20;
outsheet h Critlin Critexp Critrec using junkcrit.csv if h<=20 , comma replace ;

list h multlin1 multexp1 multexpb multrec1 multrecb if h<=20;
outsheet h multlin1 multexp1 multexpb multrec1 multrecb using junk1step.csv if h<=20, comma replace ;


list h multlin1 seylin multexp1 seyexp multrec1 seyrec ptestdiff arp if h<=20;
outsheet h multlin1 multexp1 multrec1 using junkmultse.csv if h<=20, comma replace ;


capture log close;
