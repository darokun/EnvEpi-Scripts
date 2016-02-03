*********** SAS-program: plots;

/* Forest plots for best lag data of PM2.5, PM10 and PMcoarse */

/*Important: Study names for fixed and random effects pooling have to be "pooled (fixed)"
and "pooled (random)"! They have to be at the very last and last line of the dataset!*/

/* For our exercise we use the DerSimonian-Laird (MOM)-method for the random pooling */

/* PM2.5 */
data in_bestlag1;
input OR CI_lower CI_upper IQR Study $30.;  **** Put data for 5 variables into dataset, variable study is character (author), others are numerical;
cards;
   1.0244 1.0054 1.0465 4.26 Dennekamp_2010
   1.0370 1.0070 1.0680 6.00 Ensor_2013
   1.0100 0.9000 1.1200 13.0 Levy_2001
   1.1200 1.0100 1.2500 10.0 Rosenthal_2008
   1.0700 1.0000 1.1500 7.70 Rosenthal_2013
   1.0600 1.0200 1.1000 10.0 Silverman_2010
   0.9400 0.8800 1.0200 13.8 Sullivan_2003
   1.0520 1.0000 1.0950 5.00 Wichmann_2013
   1.04984969 1.02836 1.07179 10 pooled (fixed)
   1.05054905 1.01696 1.08525 10 pooled (random)
;
run;

/* Re-calculate study estimates from IQR-increase in 10µg/m³ PM2.5-increase */

data bestlag1;
   set in_bestlag1;
   beta=log(OR)/IQR;
   lowbeta=log(CI_lower)/IQR;
   upbeta=log(CI_upper)/IQR;
   stderr=(upbeta-beta)/1.96;   **** Valid for 95%-CI;
   var=stderr**2;
   OR_new=exp(beta*10);
   CI_lower_new=exp(lowbeta*10);
   CI_upper_new=exp(upbeta*10);
run;


/* PM10 */
data in_bestlag2;
input OR CI_lower CI_upper IQR Study $30.;  **** Put data for 5 variables into dataset, variable study is character (author), others are numerical;
cards;
   0.9709 0.9485 0.9938 9.660 Dennekamp_2010
   0.8680 0.7440 1.0120 19.30 Levy_2001
   1.0600 0.9900 1.1300 14.00 Rosenthal_2013
   1.0600 0.8700 1.2700 16.51 Sullivan_2003
   1.0470 1.0070 1.0880 5.000 Wichmann_2013
   0.99021736 0.97101 1.00980 10 pooled (fixed)
   1.00995852 0.95731 1.06550 10 pooled (random)
;
run;

/* Re-calculate study estimates from IQR-increase in 10µg/m³ PM10-increase */

data bestlag2;
   set in_bestlag2;
   beta=log(OR)/IQR;
   lowbeta=log(CI_lower)/IQR;
   upbeta=log(CI_upper)/IQR;
   stderr=(upbeta-beta)/1.96;   **** Valid for 95%-CI;
   var=stderr**2;
   OR_new=exp(beta*10);
   CI_lower_new=exp(lowbeta*10);
   CI_upper_new=exp(upbeta*10);
run;


/* PMcoarse */
data in_bestlag3;
input OR CI_lower CI_upper IQR Study $30.;  **** Put data for 5 variables into dataset, variable study is character (author), others are numerical;
cards;
   0.9653 0.9226 0.9885 6.07 Dennekamp_2010
   1.0200 0.9600 1.0900 7.00 Rosenthal_2013
   1.0350 0.9970 1.0740 4.00 Wichmann_2013
   0.97225312 0.94005 1.00556 10 pooled (fixed)
   1.01125764 0.92101 1.11035 10 pooled (random)
;
run;

/* Re-calculate study estimates from IQR-increase in 10µg/m³ PMcoarse-increase */

data bestlag3;
   set in_bestlag3;
   beta=log(OR)/IQR;
   lowbeta=log(CI_lower)/IQR;
   upbeta=log(CI_upper)/IQR;
   stderr=(upbeta-beta)/1.96;   **** Valid for 95%-CI;
   var=stderr**2;
   OR_new=exp(beta*10);
   CI_lower_new=exp(lowbeta*10);
   CI_upper_new=exp(upbeta*10);
run;




%macro forest_plot (data=);
goptions reset=all;

data forest_dat;
set &data;
n=_n_;
symbol=1;
if study in ("pooled (fixed)","pooled (random)") then n=n+2;
if study="pooled (fixed)" then symbol=2;
if study="pooled (random)" then symbol=3;
run;

proc means data=forest_dat noprint;
var n;
output out=max max=max;
run;

data max;
set max;
max=max+1;
call symput('max',max);
run;


data forest_dat;
set forest_dat;
x=or_new;output;
x=ci_lower_new;symbol=4;output;
x=ci_upper_new;symbol=4;output;
drop or_new;
run;

data forest_dat;
set forest_dat;
rename x=OR_new;
run;

title angle=90 "Forest Plot";
symbol1 c=white r=4; 
symbol2 c=black i=none v=W font=marker pointlabel=("#study" h=0.7 angle=90); 
symbol3 c=red   i=none v=U font=marker pointlabel=("#study" h=0.7 angle=90);
symbol4 c=blue  i=none v=C font=marker pointlabel=("#study" h=0.7 angle=90);  
symbol5 c=black i=HILOT v=none;

axis1 label=none major=none minor=none value=none order=(0 to &max);
axis2 label=none major=none minor=none value=none;
axis3 label=(angle=90) value=(angle=90) ;

proc gplot data=forest_dat;
plot OR_new*n=symbol /vref=1 nolegend haxis=axis1 vaxis=axis2;
plot2 OR_new*n=symbol /nolegend vaxis=axis3;
run;quit;

proc datasets library=work nolist;
delete forest_dat max;
run;quit;

title;
%mend forest_plot;


ods html body="forestplotbestlag1.htm";
%forest_plot(data=bestlag1);
ods html close;

ods html body="forestplotbestlag2.htm";
%forest_plot(data=bestlag2);
ods html close;

ods html body="forestplotbestlag3.htm";
%forest_plot(data=bestlag3);
ods html close;





/* Funnel plot for best lag data of PM2.5, PM10 and PMcoarse */

/* PM2.5 */
data funnel_bestlag1;
input beta se study $30.;  **** Put data for 3 variables into dataset, variable study is character (author), others are numerical;
cards;
0.0056589379 0.0025563126 Dennekamp_2010
0.0060553215 0.0025047459 Ensor_2013
0.0007654101 0.0040572353 Levy_2001
0.0113328685 0.0056027993 Rosenthal_2008
0.0087868375 0.0047775837 Rosenthal_2013
0.0058268908 0.0018898608 Silverman_2010
-0.004483725 0.0030197438 Sullivan_2003
0.0101386229 0.0040878825 Wichmann_2013
;
run;

/* PM10 */
data funnel_bestlag2;
input beta se study $30.;  **** Put data for 3 variables into dataset, variable study is character (author), others are numerical;
cards;
-0.003057122 0.0012312768 Dennekamp_2010
-0.0073349 0.0040576328 Levy_2001
0.0041620649 0.0023304929 Rosenthal_2013
0.00352931 0.0055856065 Sullivan_2003
0.0091857864 0.0039196139 Wichmann_2013
;
run;

/* PMcoarse */
data funnel_bestlag3;
input beta se study $30.;  **** Put data for 3 variables into dataset, variable study is character (author), others are numerical;
cards;
-0.005818179 0.0019962435 Dennekamp_2010
0.0028289468 0.004837833 Rosenthal_2013
0.0086003567 0.0047179298 Wichmann_2013
;
run;


%macro funnel_plot (data=);
goptions reset=all;

data funnel_plot;
set &data;
seinv=1/se;
symbol=1;
run;

title "Funnel Plot";
axis1 label=("1/se(log(OR))");
axis2 label=("log(OR)");
symbol1  c=black i=none v=W font=marker pointlabel=("#study" h=0.7); 
proc gplot data=funnel_plot;
plot seinv*beta=symbol/href=0 vaxis=axis1 haxis=axis2 nolegend;
run;
quit;

proc datasets library=work nolist;
delete funnel_plot;
run;quit;

title;
%mend funnel_plot;


ods html body="funnelplotbestlag1.htm";
%funnel_plot(data=funnel_bestlag1);
ods html close;

ods html body="funnelplotbestlag2.htm";
%funnel_plot(data=funnel_bestlag2);
ods html close;

ods html body="funnelplotbestlag3.htm";
%funnel_plot(data=funnel_bestlag3);
ods html close;

