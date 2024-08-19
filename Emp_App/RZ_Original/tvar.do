**** TVAR.DO

*** November 21, 2016
*** Check for explosive roots in TVAR and estimate

*** Uses Gordon-Krenn (GK) transformation (i.e. variables are divided by potential)

***
*** Requires:
***     rzdat.xlsx, updated November 21, 2016
********************************************************************************

 #delimit;

drop _all;
clear all;

set more 1;
set matsize 800;

capture log close;
log using tvar_results.log, replace;

*******************************************************************************;
** RAW DATA IMPORTATION AND DATA SETUP;
*******************************************************************************;

*import excel rzdat.xlsx, sheet("rzdat") firstrow;

insheet using rzdatnew.csv;

drop if quarter<1889.75;

gen qdate = q(1889q4) + _n-1;
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


* Calculate mean, median, etc. of slack duration;
/*
tsspell slack, end(hrend) seq(hrseq);
summ hrseq if hrend & slack==1, detail;
tab hrseq if hrend & slack==1;
*/
*******************************************************************************;
/* Options for VARs below

change sample according to:

no "if" is full sample
if L.zlb ==1 or L.zlb==0
if L.slack ==1 or L.slack==0;

add taxes:  var newsy g y taxy or var g y taxy;

add trends:  var newsy g y taxy if ..., lags(1/4) exog(t t2);

*******************************************************************************/

/* Military News identification;*/

var newsy g  y if L.zlb== 1, lags(1/4);
varstable, graph name(stability);
irf create irf, step(20) set(irf, replace);
irf table oirf, impulse(newsy) response(newsy g y) ;
irf graph oirf, impulse(newsy) response(newsy g  y) byopts(rescale) saving(newszlb.gph, replace);

var newsy g  y if L.zlb== 0, lags(1/4);
varstable, graph name(stability2);
irf create irf, step(20) set(irf, replace);
irf table oirf, impulse(newsy) response(newsy g y) ;
irf graph oirf, impulse(newsy) response(newsy g  y) byopts(rescale) saving(newsnozlb.gph, replace);


var newsy g  y if L.recession== 1, lags(1/4);
varstable, graph name(stability3);
irf create irf, step(20) set(irf, replace);
irf table oirf, impulse(newsy) response(newsy g y) ;
irf graph oirf, impulse(newsy) response(newsy g  y) byopts(rescale) saving(newsslack.gph, replace);


var newsy g  y if L.recession== 0, lags(1/4);
varstable, graph name(stability4);
irf create irf, step(20) set(irf, replace);
irf table oirf, impulse(newsy) response(newsy g y) ;
irf graph oirf, impulse(newsy) response(newsy g  y) byopts(rescale) saving(newsnoslack.gph, replace);


* Blanchard-Perotti identification;


var g y if L.zlb == 1 , lags(1/4); 
varstable, graph name(stability5);
irf create irfbpzlb, step(20) set(bpirf, replace);
irf table oirf, impulse(g) response( g y) ;
irf graph oirf, impulse(g) response(g y) byopts(rescale) saving(bpzlb.gph, replace);


var g y if L.zlb == 0, lags(1/4) ; 
varstable, graph name(stability6);
irf create irfbpnozlb, step(20) set(bpirf, replace);
irf table oirf, impulse(g) response( g y) ;
irf graph oirf, impulse(g) response(g y) byopts(rescale) saving(bpnozlb.gph, replace);


var g y if L.recession == 1, lags(1/4) ; 
irf create irfbpslack, step(20) set(bpirf, replace);
irf table oirf, impulse(g) response( g y) ;
irf graph oirf, impulse(g) response(g y) byopts(rescale) saving(bpslack.gph, replace);

var g y if L.recession == 0, lags(1/4) ; 
irf create irfbpnoslack, step(20) set(bpirf, replace);
irf table oirf, impulse(g) response( g y) ;
irf graph oirf, impulse(g) response(g y) byopts(rescale) saving(bpnoslack.gph, replace);


/*
graph combine bplin.gph bpzlb.gph bpnozlb.gph bpslack.gph bpnoslack.gph;
*/



capture log close;

 
