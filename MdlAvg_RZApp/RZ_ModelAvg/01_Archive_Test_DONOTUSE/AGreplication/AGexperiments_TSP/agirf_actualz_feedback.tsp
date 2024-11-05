name ramey "Auerbach-Gorodnichenk Model and Data: This program allows actual state dynamics and partial feedback";

freq q; smpl 1947:1 2009:3;

??    READ IN AUERBACH-GORODNICHENKO DATA

read(file="agtspdat.xls");


??    CREATE VARIABLES

trend t; t2 = t**2;


gen ly = lrgdp;
gen lg = lrgovag;
gen lt = lrnettax;

gen lym = 5.6146*ly; ? AG'S AD HOC CONVERSION FACTOR

gen lgap = ly - log(rgdp_pot);  ? CREATE LOG OUTPUT GAP


? CREATE AG'S ECONOMIC STATE INDICATOR

smpl 1948:1 2008:4;

z = 100*(ly(3) - ly(-4))/7;  ? AG'S CENTERED MOVING AVERAGE OF GROWTH

?The mean of z is approx. 0.8 and std is 0.5.  AG specically use those numbers so we do too

znorm = (z - 0.8)/0.5;

fznorm = exp(-1.5*znorm)/(1 + exp(-1.5*znorm));

gen fznormlag = fznorm(-1);


? SIMULATION SET UP

smpl 1947:1 2009:4;


??     SET THE RESIDUALS TO 0 TO START WITH

gen gres = 0;
gen tres = 0;
gen yres = 0;

gshock = 0;


??    CREATE THE EQUATIONS FROM THE STVAR MODEL

frml geq lg = g0 +
  (1-fznorm(-1))*[ ge1*lg(-1) + ge2*lt(-1) + ge3*ly(-1)
                 + ge4*lg(-2) + ge5*lt(-2) + ge6*ly(-2)
                 + ge7*lg(-3) + ge8*lt(-3) + ge9*ly(-3)]
  +   fznorm(-1)*[ gr1*lg(-1) + gr2*lt(-1) + gr3*ly(-1)
                 + gr4*lg(-2) + gr5*lt(-2) + gr6*ly(-2)
                 + gr7*lg(-3) + gr8*lt(-3) + gr9*ly(-3)]
  + gz1*z(-1) + gz2*z(-2) + gz3*z(-3) + gz4*z(-4)
  + gg*(gres + gshock);


frml teq lt = t0 +
  (1-fznorm(-1))*[ te0*(gres + gshock) + te1*lg(-1) + te2*lt(-1) + te3*ly(-1)
                 + te4*lg(-2) + te5*lt(-2) + te6*ly(-2)
                 + te7*lg(-3) + te8*lt(-3) + te9*ly(-3)]
  +   fznorm(-1)*[ tr0*(gres + gshock) + tr1*lg(-1) + tr2*lt(-1) + tr3*ly(-1)
                 + tr4*lg(-2) + tr5*lt(-2) + tr6*ly(-2)
                 + tr7*lg(-3) + tr8*lt(-3) + tr9*ly(-3)]
  + tz1*z(-1) + tz2*z(-2) + tz3*z(-3) + tz4*z(-4)
  + tt*tres;


frml yeq ly = y0 +
  (1-fznorm(-1))*[ ye0*(gres + gshock) + ye1*lg(-1) + ye2*lt(-1) + ye3*ly(-1)
                 + ye4*lg(-2) + ye5*lt(-2) + ye6*ly(-2)
                 + ye7*lg(-3) + ye8*lt(-3) + ye9*ly(-3)]
  +   fznorm(-1)*[ yr0*(gres + gshock) + yr1*lg(-1) + yr2*lt(-1) + yr3*ly(-1)
                 + yr4*lg(-2) + yr5*lt(-2) + yr6*ly(-2)
                 + yr7*lg(-3) + yr8*lt(-3) + yr9*ly(-3)]
   + yz1*z(-1) + yz2*z(-2) + yz3*z(-3) + yz4*z(-4)
  + yy*yres;



param g0,0,t0,0,y0,0;

const gg,0,tt,0,yy,0;


??  SET THE PARAMETERS EQUAL TO THOSE WE ESTIMATED IN MATLAB USING AG PROGRAMS AND DATA

? Updated 11:47am June 9 2014

const gr1,0.83279,tr1,0.061513,yr1,0.0067301
gr2,-0.029452,tr2,0.74741,yr2,-0.0039289
gr3,0.10811,tr3,0.96164,yr3,1.0791
gr4,0.14329,tr4,0.50262,yr4,0.091612
gr5,-0.071073,tr5,0.53472,yr5,0.097508
gr6,0.12981,tr6,-1.24,yr6,-0.237
gr7,-0.027595,tr7,-0.37823,yr7,-0.0782
gr8,0.13889,tr8,-0.4037,yr8,-0.093321
gr9,-0.22673,tr9,0.24332,yr9,0.13929
ge1,1.5062,te1,-0.36795,ye1,-0.036976
ge2,0.018932,te2,0.90514,ye2,0.039461
ge3,-0.11219,te3,1.1282,ye3,0.85782
ge4,-0.37924,te4,-0.13209,ye4,0.048564
ge5,0.26683,te5,-0.34983,ye5,-0.13483
ge6,-0.46502,te6,-1.614,ye6,-0.15617
ge7,-0.15577,te7,0.38383,ye7,-0.028043
ge8,-0.26841,te8,0.40337,ye8,0.10029
ge9,0.58067,te9,0.62628,ye9,0.30582
gz1,0.039832,tz1,0.0324,yz1,0.011526
gz2,-0.007058,tz2,-0.030025,yz2,-0.0093221
gz3,0.0062417,tz3,0.063996,yz3,0.036003
gz4,-0.00027767,tz4,-0.03079,yz4,-0.030354;


const te0,0.15816,ye0,0.10008,tr0,0.26247,yr0,0.095841;  ?FROM MATLAB IMPULSE RESPONSES 


smpl 1949:3 2008:4;

const g0,-0.0112, t0,-0.1282, y0, 0.0207;  ? FROM MATLAB ESTIMATION


?  THE FOLLOWING LINES CALCULATE THE RESIDUALS AND SET THE ONES IN THE EQUATIONS EQUAL TO THEM

smpl 1949:3 2008:4;

gen geq gfit; gen gres = lg - gfit;
genr teq tfit; gen tres = lt - tfit;
genr yeq yfit; gen yres = ly - yfit;


smpl 1949:3 2008:4;

msd gres tres yres;  sdgres = @stddev(1);

ident yieq lym = 5.6146*ly;


?SET COEFFICIENTS ON OWN SHOCKS TO 1

const gg,1, tt,1,yy,1;



? TEST THAT WE ARE USING THE RIGHT RESIDUALS BY COMPARING SIMULATED TO ACTUAL;


smpl 1949:3 2008:4;

 siml(dynam,noprnsim,endog=(lg lt ly lym ),tag=ss) geq teq yeq yieq ;

title "CHECK THAT THE SIMULATED IS IDENTICAL TO THE ACTUAL";

print lg lgss ly lyss;  

smpl 1947:1 2008:4;

? -A- SUFFIX STANDS FOR ACTUAL

dot ly lg lt lym z znorm fznorm;
  gen .a = .;
enddot;


smpl 1949:3 2008:4;

trend t; 

?NOW CONSIDER SHOCKS;


??  THIS PART CALCULATES THE IMPULSE RESPONSES ALLOWING FZNORM TO VARY AS IT DOES IN THE DATA
??   WE COMPARE THE PATH WITH THE SHOCK TO THE ONE THAT THE ECONOMY WAS ON AND CALCULATE 
??     THE MULTIPLIERS FROM THAT.

? WE CALCULATE MULTIPLIERS FOR 8 QUARTERS AND 20 QUARTERS


? gshock = 0.01;


ident zeq z = 100*(ly(3) - ly(-4))/7;

ident znormeq znorm = (z - 0.8)/0.5;

ident fzeq fznorm = exp(-1.5*znorm)/(1 + exp(-1.5*znorm));


do i = 1 to 238;

   smpl 1949:3 2008:4;
   gshock = 0;
   select t = i; gshock = sdgres; 
   select t>=i & t<i + 20;

   siml(dynam,silent,endog=(lg lt ly lym z znorm fznorm),tag=s)
         geq teq yeq yieq zeq znormeq fzeq;

?   print fznorm fznorma fznorms;
?   print ly lya lys;
?   print lg lga lgs;

?  dot ly lg;
?    plot(noorigin) .a,a,.s,s;
?  enddot;

  dlyms = lyms - lyma;
  dlgs = lgs - lga;
  msd(silent) dlyms dlgs;
  set m20 = @mean(1)/@mean(2);

  select t>=i & t<i + 8;
  msd(silent) dlyms dlgs;
  set m8 = @mean(1)/@mean(2);


  select t = i;  tindex = i; mult8 = m8; mult20 = m20;

enddo;

smpl 2004:2 2008:4; mult20 = @miss;
smpl 2007:2 2008:4; mult8 = @miss;

smpl 1949:3 2008:4;


title "MULTIPLIERS FOR AG EXTREME RECESSIONS";

select fznorm(-1)>=0.95; msd mult8 mult20 fznorm(-1);

title "MULTIPLIERS FOR AG EXTREME BOOMS";

select fznorm(-1)<=0.05; msd mult8 mult20 fznorm(-1);



stop; end;


end;
