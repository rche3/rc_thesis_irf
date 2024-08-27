/* empirical_longhor.prg
See "explain_longhor.pdf".

Let     Ds(t) = pct change in quarterly exchange rate
        yy(t) = (Ds(t)+DS(t+1)+...+DS(t+q)/(q+1) = average value of q+1 period change in exchange rate from t-1 to t+q
        x(t-1) = q+1 period nominal interest differential

Get bias adjusted estimate of beta in the interest parity regression

        yy(t) = const. + beta*x(t-1) + eta(t)

Data are in "empirical_example.xls".  For the UK
        Ds is called uk_er_3m
        yy is called  uk_er_6m for q=1, or uk_er_1y for q=4,
        x is called uk_id_3m, uk_id_6m, etc.
Names for Japan and Canada are the same, with Canada or jap replacing uk.

Note the odd dating of yy(t) in the spreadsheet: the entry for, say, 1979:2 and q=1 is [Ds(79:2)+Ds(79:3)]/2
and thus is realized in 79:3 even though it appears in the 79:2 row in the spreadsheet.  For example,
for the UK and q=1, uk_er_6m in 1979:2 is 12.08 = (18.81+5.35)/2, where 18.81 and 5.35 are the 79:2 and 79:3
values of uk_er_3m.

Needs to be run separately for each currency and horizon in the data set.

Before running, change "open data" to open the desired ".inp" file.  The .inp filenames specifies
        *currency and horizon
        *whether univariate AR in x(nvar=1) or bivariate VAR in x and Ds (nvar=2) is used to compute the relevant moments
        *number of lags to use in the AR or VAR (nlag, only allowed values are nlag=1 or nlag=2)
The names of the .inp file have the structure
        <country>var<q>.inp

The .inp files supplied VAR with two lags.

For example Canadavar3.inp is the following:
 1979     1      2011  2     4    2       2      Canada_id_1y   Canada_er_1y  Canada_er_3m          Canadavar3
start_y start_q end_y end_q  nq+1  nlag    nvar    i_series      er_series     er_one_period     output file name

And UKvar19.inp
 1979     1      2011  3     20    2       2      UK_id_5y      UK_er_5y     UK_er_3m          UKvar19
start_y start_q end_y end_q  nq+1  nlag    nvar    i_series      er_series    er_one_period     output file name

As noted above, in the spreadsheet yy(t) is realized in period t+q.  To invoke longhor, the code sets

        yseries(t+q) = yy(t)

so that yseries(t) is realized in period t.
*/

/*
Structure of the program

1. "source" to include routines invoked by the program, along with a procedure to compute the s.e.
2a. declare some variables
2b. read in data
3a. declare variables needed by longhor
3b. invoke longhor
4. compute s.e.
5. display results

*/

* 1. called routines
source(noecho) c:\a\myfiles\bias\replication_files\proc_vbias.src
source(noecho) c:\a\myfiles\bias\replication_files\proc_vb_ma0.src
source(noecho) c:\a\myfiles\bias\replication_files\proc_vb_maq.src
source(noecho) c:\a\myfiles\bias\replication_files\longhor.src
source(noecho) c:\a\myfiles\bias\replication_files\longhor1.src
source(echo) c:\a\myfiles\bias\replication_files\nwbandwidth.src

*
Procedure estimate_IP yseries xseries start_IP end_IP qnow betahat_ip se_ip

  type series yseries xseries
  type int start_IP end_IP qnow
  type real *betahat_ip *se_ip
* parameters for bandwidth calculation in newey-west
  declare real s0 s1
  declare int m
  *eta_ols_hat --- least square residuals
  declare series[real] eta_ols_hat xnw

  linreg(noprint) yseries start_IP end_IP eta_ols_hat
  # constant xseries{qnow+1}

  set xnw start_IP-1 end_IP-1 = xseries(t-qnow)

  * compute bandwidth m for the IP regression
  @nwbandwidth xnw eta_ols_hat start_IP end_IP qnow s0 s1 m
  display "estimate_IP: s0, s1, m:" s0 s1 m

  linreg(noprint,robusterrors,lwindow=neweywest,lags=m) yseries start_IP end_IP eta
  # constant xseries{qnow+1}

  compute betahat_ip = %BETA(2)
  compute se_ip = %STDERRS(2)

end estimate_IP


* 2a. declare variables and read in control parameters
calendar(q) 1970:1
allocate 2011:4

declare int narlag nvar first last nqplus1 nq
declare int start_y start_q end_y end_q
declare int start_IP end_IP qnow
declare real betahat_ip se_ip
declare label i_series er_series er_one_period out_file
*************
*** change this to select the currency and horizon
open data UKvar39.inp
***
*************
read(unit=data) start_y start_q end_y end_q nqplus1 narlag nvar i_series er_series er_one_period out_file
close data

* 2b. read in data

declare series[real] x y y_lag v

comp temp_first=%cal(start_y,start_q)
comp temp_last=%cal(end_y,end_q)
disp %datelabel(temp_first) %datelabel(temp_last) nqplus1 narlag i_series er_series er_one_period

open data "example_data.xlsx"
data(org=columns,format=xlsx,sheet="Sheet1",left=2) 1970:1 2011:4 %s(i_series) %s(er_series) %s(er_one_period)
close data
compute nq=nqplus1-1
set x temp_first temp_last = %s(i_series)
set yy temp_first temp_last = %s(er_series)
set Ds temp_first temp_last = %s(er_one_period)
print(dates) 1979:1 1979:4 x yy Ds

* 3a. The following variables are for longhor

declare int nk nWseries
compute nk=1, nWseries=1
declare vector[real] vbias(nk) betahat(nk) betahat_adj(nk)
declare vector[series] Wseries(nWseries)
comp first=temp_first+1+nq
comp last=temp_last
disp %datelabel(first) %datelabel(last) nqplus1 narlag i_series er_series er_one_period

set yseries / = yy(t-nq)
set xseries / = x(t)
set Wseries(1) / = Ds(t)

* 3b. use longhor

if nvar .eq. 1
  execute longhor1 yseries xseries first last nq nk narlag vbias betahat betahat_adj
else
  execute longhor yseries xseries Wseries nWseries first last nq nk narlag vbias betahat betahat_adj
end

* 4. Get standard error
*Awkward: to get s.e., need to re-estimate IP.
execute estimate_IP yseries xseries first last nq betahat_ip se_ip
comp beta_diff = betahat(1)-betahat_ip
comp tol = 0.000001

* 5. Display results

*do is just to put us in compiler mode, so that all "display" output is together
do i=1,1
if abs(beta_diff) .gt. tol; begin
  disp "betahat estimated in code to get s.e. does not match longhor betahat!!"
  comp se_ip = %na
end
disp " "
disp " "
disp "series for q+1 period exchange rate is" er_series " with q= " nq "and nvar= " nvar
disp "b"       @20 ###.### vbias
disp "b/T"  @20 ###.### vbias/Tnow
disp "betahat"     @20 ###.### betahat
disp "betahat_adj" @20 ###.### betahat_adj
disp @20 "(" @-1 ##.### se_ip ")"
disp @5 "betahat"  @16 "b/T" @25 "betahat_adj"
disp @5 ##.## betahat @15 ##.## vbias/Tnow @25 betahat_adj
disp @4 "(" @-1 ##.## se_ip ")"

disp " "
end do

end



