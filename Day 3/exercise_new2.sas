*********** SAS-program: exercise;


/* Read estimates from selected papers into a SAS dataset */

**** Analysis of PM2.5 best lag: RRs or ORs, CIs and IQRs (from 8 papers);
data in_estimates_bestlag1; 
   input study oddsratio lowerCI upperCI iqr author $30.;   **** Put data for 6 variables into dataset, variable study counts papers (1-8), variable author is character, others are numerical;
   cards;                                            
   1 1.0244 1.0054 1.0465 4.26 Dennekamp_2010
   2 1.0370 1.0070 1.0680 6.00 Ensor_2013
   3 1.0100 0.9000 1.1200 13.0 Levy_2001
   4 1.1200 1.0100 1.2500 10.0 Rosenthal_2008
   5 1.0700 1.0000 1.1500 7.70 Rosenthal_2013
   6 1.0600 1.0200 1.1000 10.0 Silverman_2010
   7 0.9400 0.8800 1.0200 13.8 Sullivan_2003
   8 1.0520 1.0000 1.0950 5.00 Wichmann_2013
   ;                             
run; 

data estimates_bestlag1;
   set in_estimates_bestlag1;
   beta=log(oddsratio)/iqr;
   lowbeta=log(lowerCI)/iqr;
   upbeta=log(upperCI)/iqr;
   stderr=(upbeta-beta)/1.96;   **** Valid for 95%-CI;
   var=stderr**2;
run;

/* Apply macro: let SAS run the macro "pooling" once before you move on! */

**** Pooling of best lag PM2.5 estimates;
%pooling (estimates_bestlag1);
quit; 

/* Convert pooled estimate and its SE into OR and 95%-CI per PM2.5-increase of 10µg/m³ */

data all_or;
   set all;
   fe_lowci=exp(10*fe-1.96*10*fe_se);
   fe_upci=exp(10*fe+1.96*10*fe_se);
   fe=exp(10*fe);
   drop fe_se;   
   re_dl_lowci=exp(10*re_dl-1.96*10*se_re_dl);
   re_dl_upci=exp(10*re_dl+1.96*10*se_re_dl);
   re_dl=exp(10*re_dl);
   drop se_re_dl;   
   re_reml_lowci=exp(10*re_reml-1.96*10*se_re_reml);
   re_reml_upci=exp(10*re_reml+1.96*10*se_re_reml);
   re_reml=exp(10*re_reml);
   drop se_re_reml;   
   re_ml_lowci=exp(10*re_ml-1.96*10*se_re_ml);
   re_ml_upci=exp(10*re_ml+1.96*10*se_re_ml);
   re_ml=exp(10*re_ml);
   drop se_re_ml;
run;

/* Create nice output dataset for SAS and Excel */

data all_or;
   retain daten prob_hom heterogene fe fe_lowci fe_upci prob_fe tau_2 re_dl re_dl_lowci 
          re_dl_upci prob_re_dl re_reml re_reml_lowci re_reml_upci prob_re_reml re_ml re_ml_lowci 
          re_ml_upci prob_re_ml;    **** Put variables in correct order;
   set all_or;
   if prob_hom < 0.10 then heterogene=1;     **** Use p=0.10 as conservative limit for test of homogeneity;
      else heterogene=0;
   if heterogene=0 then do; pooled=fe; pooled_lowci=fe_lowci; pooled_upci=fe_upci; pooled_prob=prob_fe; end;			         **** Use fixed effects estimates for homogeneity;
   if heterogene=1 then do; pooled=re_dl; pooled_lowci=re_dl_lowci; pooled_upci=re_dl_upci; pooled_prob=prob_re_dl; end;         **** Use random effects estimates (DerSimonian-Laird) for heterogeneity;
   label daten="Data" fe="Fixed Effects Estimate" fe_lowci="Lower CI for FE Estimate" fe_upci="Upper CI for FE Estimate" 
         prob_hom="p-Value Test for Homogeneity" heterogene="Heterogeneity" 
         re_dl="MOM RE Estimate (DerSimonian-Laird)"  re_dl_lowci="Lower CI for MOM RE Estimate" re_dl_upci="Upper CI for MOM RE Estimate"
		 re_reml="REML RE Estimate (Normand)"  re_reml_lowci="Lower CI for REML RE Estimate"  re_reml_upci="Upper CI for REML RE Estimate"
		 re_ml="ML RE Estimate (van Houwelingen)" re_ml_lowci="Lower CI for ML RE Estimate" re_ml_upci="Upper CI for ML RE Estimate"
         prob_fe="p-Value for FE Estimate" prob_re_dl="p-Value for MOM RE Estimate" prob_re_ml="p-Value for ML RE Estimate" prob_re_reml="p-Value for REML RE Estimate"  
         pooled="Pooled Estimate" pooled_lowci="Lower CI for Pooled Estimate" pooled_upci="Upper CI for Pooled Estimate" pooled_prob="p-Value for Pooled Estimate" 
         tau_2="tau2";                                            
run;

/* Output of pooled best lag PM2.5 estimate */
/* Input of a path needed, e.g. 'D:\pooled_bestlag1_estimate.xls' */

ods html file='/home/msc1426/Dokumente/EnvEpi-Scripts/Day 3/pooled_bestlag1_estimate.xls' headtext="<style>   
td {mso-number-format:\@}</style>";       **** Output to excel sheet;
title1 "Pooled Best Lag Estimate PM2.5";

proc print data=all_or noobs label; 
run;

ods html close;
quit; 



/* Funnel plot regression for best lag PM2.5 estimates */

data estimates_bestlag1;
   set estimates_bestlag1;
   precision=1/stderr;
   SND=beta/stderr;
run;

ods html body="funnelplotregbestlag1.htm";
proc reg data=estimates_bestlag1;
    model SND=precision; 
    plot SND*precision; 
run;
ods html close;
quit;






**** Analysis of PM10 best lag: RRs or ORs, CIs and IQRs (from 5 papers);
data in_estimates_bestlag2; 
   input study oddsratio lowerCI upperCI iqr author $30.;   **** Put data for 6 variables into dataset, variable study counts papers (1-5), variable author is character, others are numerical;
   cards;                                            
   1 0.9709 0.9485 0.9938 9.660 Dennekamp_2010
   2 0.8680 0.7440 1.0120 19.30 Levy_2001
   3 1.0600 0.9900 1.1300 14.00 Rosenthal_2013
   4 1.0600 0.8700 1.2700 16.51 Sullivan_2003
   5 1.0470 1.0070 1.0880 5.000 Wichmann_2013
   ;                           
run; 

data estimates_bestlag2;
   set in_estimates_bestlag2;
   beta=log(oddsratio)/iqr;
   lowbeta=log(lowerCI)/iqr;
   upbeta=log(upperCI)/iqr;
   stderr=(upbeta-beta)/1.96;   **** Valid for 95%-CI;
   var=stderr**2;
run;

/* Apply macro: let SAS run the macro "pooling" once before you move on! */

**** Pooling of best lag PM10 estimates;
%pooling (estimates_bestlag2);
quit; 

/* Convert pooled estimate and its SE into OR and 95%-CI per PM10-increase of 10µg/m³ */

data all_or;
   set all;
   fe_lowci=exp(10*fe-1.96*10*fe_se);
   fe_upci=exp(10*fe+1.96*10*fe_se);
   fe=exp(10*fe);
   drop fe_se;   
   re_dl_lowci=exp(10*re_dl-1.96*10*se_re_dl);
   re_dl_upci=exp(10*re_dl+1.96*10*se_re_dl);
   re_dl=exp(10*re_dl);
   drop se_re_dl;   
   re_reml_lowci=exp(10*re_reml-1.96*10*se_re_reml);
   re_reml_upci=exp(10*re_reml+1.96*10*se_re_reml);
   re_reml=exp(10*re_reml);
   drop se_re_reml;   
   re_ml_lowci=exp(10*re_ml-1.96*10*se_re_ml);
   re_ml_upci=exp(10*re_ml+1.96*10*se_re_ml);
   re_ml=exp(10*re_ml);
   drop se_re_ml;
run;

/* Create nice output dataset for SAS and Excel */

data all_or;
   retain daten prob_hom heterogene fe fe_lowci fe_upci prob_fe tau_2 re_dl re_dl_lowci 
          re_dl_upci prob_re_dl re_reml re_reml_lowci re_reml_upci prob_re_reml re_ml re_ml_lowci 
          re_ml_upci prob_re_ml;    **** Put variables in correct order;
   set all_or;
   if prob_hom < 0.10 then heterogene=1;     **** Use p=0.10 as conservative limit for test of homogeneity;
      else heterogene=0;
   if heterogene=0 then do; pooled=fe; pooled_lowci=fe_lowci; pooled_upci=fe_upci; pooled_prob=prob_fe; end;			         **** Use fixed effects estimates for homogeneity;
   if heterogene=1 then do; pooled=re_dl; pooled_lowci=re_dl_lowci; pooled_upci=re_dl_upci; pooled_prob=prob_re_dl; end;         **** Use random effects estimates (DerSimonian-Laird) for heterogeneity;
   label daten="Data" fe="Fixed Effects Estimate" fe_lowci="Lower CI for FE Estimate" fe_upci="Upper CI for FE Estimate" 
         prob_hom="p-Value Test for Homogeneity" heterogene="Heterogeneity" 
         re_dl="MOM RE Estimate (DerSimonian-Laird)"  re_dl_lowci="Lower CI for MOM RE Estimate" re_dl_upci="Upper CI for MOM RE Estimate"
		 re_reml="REML RE Estimate (Normand)"  re_reml_lowci="Lower CI for REML RE Estimate"  re_reml_upci="Upper CI for REML RE Estimate"
		 re_ml="ML RE Estimate (van Houwelingen)" re_ml_lowci="Lower CI for ML RE Estimate" re_ml_upci="Upper CI for ML RE Estimate"
         prob_fe="p-Value for FE Estimate" prob_re_dl="p-Value for MOM RE Estimate" prob_re_ml="p-Value for ML RE Estimate" prob_re_reml="p-Value for REML RE Estimate"  
         pooled="Pooled Estimate" pooled_lowci="Lower CI for Pooled Estimate" pooled_upci="Upper CI for Pooled Estimate" pooled_prob="p-Value for Pooled Estimate" 
         tau_2="tau2";                                            
run;

/* Output of pooled best lag PM10 estimate */
/* Input of a path needed, e.g. 'D:\pooled_bestlag2_estimate.xls' */

ods html file='path\pooled_bestlag2_estimate.xls' headtext="<style>   
td {mso-number-format:\@}</style>";       **** Output to excel sheet;
title1 "Pooled Best Lag Estimate PM10";

proc print data=all_or noobs label; 
run;

ods html close;
quit;



/* Funnel plot regression for best lag PM10 estimates */

data estimates_bestlag2;
   set estimates_bestlag2;
   precision=1/stderr;
   SND=beta/stderr;
run;

ods html body="funnelplotregbestlag2.htm";
proc reg data=estimates_bestlag2;
    model SND=precision; 
    plot SND*precision; 
run;
ods html close;
quit;





**** Analysis of PMcoarse best lag: RRs or ORs, CIs and IQRs (from 3 papers);
data in_estimates_bestlag3; 
   input study oddsratio lowerCI upperCI iqr author $30.;   **** Put data for 6 variables into dataset, variable study counts papers (1-3), variable author is character, others are numerical;
   cards;                                            
   1 0.9653 0.9226 0.9885 6.07 Dennekamp_2010
   2 1.0200 0.9600 1.0900 7.00 Rosenthal_2013
   3 1.0350 0.9970 1.0740 4.00 Wichmann_2013
   ;                           
run; 

data estimates_bestlag3;
   set in_estimates_bestlag3;
   beta=log(oddsratio)/iqr;
   lowbeta=log(lowerCI)/iqr;
   upbeta=log(upperCI)/iqr;
   stderr=(upbeta-beta)/1.96;   **** Valid for 95%-CI;
   var=stderr**2;
run;

/* Apply macro: let SAS run the macro "pooling" once before you move on! */

**** Pooling of best lag PMcoarse estimates;
%pooling (estimates_bestlag3);
quit; 

/* Convert pooled estimate and its SE into OR and 95%-CI per PMcoarse-increase of 10µg/m³ */

data all_or;
   set all;
   fe_lowci=exp(10*fe-1.96*10*fe_se);
   fe_upci=exp(10*fe+1.96*10*fe_se);
   fe=exp(10*fe);
   drop fe_se;   
   re_dl_lowci=exp(10*re_dl-1.96*10*se_re_dl);
   re_dl_upci=exp(10*re_dl+1.96*10*se_re_dl);
   re_dl=exp(10*re_dl);
   drop se_re_dl;   
   re_reml_lowci=exp(10*re_reml-1.96*10*se_re_reml);
   re_reml_upci=exp(10*re_reml+1.96*10*se_re_reml);
   re_reml=exp(10*re_reml);
   drop se_re_reml;   
   re_ml_lowci=exp(10*re_ml-1.96*10*se_re_ml);
   re_ml_upci=exp(10*re_ml+1.96*10*se_re_ml);
   re_ml=exp(10*re_ml);
   drop se_re_ml;
run;

/* Create nice output dataset for SAS and Excel */

data all_or;
   retain daten prob_hom heterogene fe fe_lowci fe_upci prob_fe tau_2 re_dl re_dl_lowci 
          re_dl_upci prob_re_dl re_reml re_reml_lowci re_reml_upci prob_re_reml re_ml re_ml_lowci 
          re_ml_upci prob_re_ml;    **** Put variables in correct order;
   set all_or;
   if prob_hom < 0.10 then heterogene=1;     **** Use p=0.10 as conservative limit for test of homogeneity;
      else heterogene=0;
   if heterogene=0 then do; pooled=fe; pooled_lowci=fe_lowci; pooled_upci=fe_upci; pooled_prob=prob_fe; end;			         **** Use fixed effects estimates for homogeneity;
   if heterogene=1 then do; pooled=re_dl; pooled_lowci=re_dl_lowci; pooled_upci=re_dl_upci; pooled_prob=prob_re_dl; end;         **** Use random effects estimates (DerSimonian-Laird) for heterogeneity;
   label daten="Data" fe="Fixed Effects Estimate" fe_lowci="Lower CI for FE Estimate" fe_upci="Upper CI for FE Estimate" 
         prob_hom="p-Value Test for Homogeneity" heterogene="Heterogeneity" 
         re_dl="MOM RE Estimate (DerSimonian-Laird)"  re_dl_lowci="Lower CI for MOM RE Estimate" re_dl_upci="Upper CI for MOM RE Estimate"
		 re_reml="REML RE Estimate (Normand)"  re_reml_lowci="Lower CI for REML RE Estimate"  re_reml_upci="Upper CI for REML RE Estimate"
		 re_ml="ML RE Estimate (van Houwelingen)" re_ml_lowci="Lower CI for ML RE Estimate" re_ml_upci="Upper CI for ML RE Estimate"
         prob_fe="p-Value for FE Estimate" prob_re_dl="p-Value for MOM RE Estimate" prob_re_ml="p-Value for ML RE Estimate" prob_re_reml="p-Value for REML RE Estimate"  
         pooled="Pooled Estimate" pooled_lowci="Lower CI for Pooled Estimate" pooled_upci="Upper CI for Pooled Estimate" pooled_prob="p-Value for Pooled Estimate" 
         tau_2="tau2";                                            
run;

/* Output of pooled best lag PMcoarse estimate */
/* Input of a path needed, e.g. 'D:\pooled_bestlag3_estimate.xls' */

ods html file='/home/msc1426/Dokumente/EnvEpi-Scripts/Day 3/pooled_bestlag3_estimate.xls' headtext="<style>   
td {mso-number-format:\@}</style>";       **** Output to excel sheet;
title1 "Pooled Best Lag Estimate PMcoarse";

proc print data=all_or noobs label; 
run;

ods html close;
quit;



/* Funnel plot regression for best lag PMcoarse estimates */

data estimates_bestlag3;
   set estimates_bestlag3;
   precision=1/stderr;
   SND=beta/stderr;
run;

ods html body="funnelplotregbestlag3.htm";
proc reg data=estimates_bestlag3;
    model SND=precision; 
    plot SND*precision; 
run;
ods html close;
quit;



