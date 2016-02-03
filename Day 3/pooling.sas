*********** SAS-program: pooling;

/* Purpose: calculate pooled estimates using fixed or random effects methods and test of homogeneity */
   
/* Basis:   Four Papers (first two are "Tutorials in Biostatistics" in Statistics in Medicine):

			    SLT Normand "Meta-Analysis: Formulating, Evaluating, Combining, and Reporting",
			    Statistics in Medicine 18, 321-359 (1999)

			    HC van Houwelingen, LR Arends, T Stijnen "Advanced methods in meta-analysis: 
			    multivariate approach and meta-regression",
			    Statistics in Medicine 21, 589-624 (2002)

                R DerSimonian, N Laird "Meta-analysis in clinical trials", 
                Controlled Clinical Trials 7, 177-188 (1986)

                R DerSimonian, R Kacker "Random-effects model for meta-analysis of clinical 
                trials: an update",
                Contemporary Clinical Trials 28, 105-114 (2007)*/


/*    A fixed-effects model assumes that each of the k study summary statistics (Yi) is a 
      realization from a population of study estimates with common mean theta (T):
      E(Yi)=T and si**2=var(Yi), Yi from N(T, si**2). 

      The random-effects framework postulates that each study summary statistic (Yi) is a draw 
      from a distribution with a study-specific mean Ti and variance si**2, Yi from N(Ti, si**2)
      and each study-specific mean Ti is assumed to be a draw from some superpopulation of 
      effects with mean T and variance tau**2: Ti from N(T, tau**2). 
      T and tau**2 are hyperparameters and represent the average effect and inter-study variation.
      After averaging over the study-specific effects, the distribution of each study summary Yi
      is normal with mean T and varicane si**2+tau**2: Yi from N(T, si**2+tau**2). As in the fixed-effects 
      model the parameter of central interest is T; however, the between-study variation tau**2 
      plays an important role and must also be estimated. It is also possible to derive estimates
      for the study-specific effects Ti.    

      There are at least three sources of variation to consider before combining summary 
      statistics across studies:
      1. Sampling error (e.g. varying sample sizes resulting in study summaries estimated with 
         varying degrees of freedom) - solution: weighted averages
      2. Study-level characteristics may differ between studies (e.g. studies in different types
         of settings causing different effects) - solution: stratification before combining estimates
      3. Inter-study variation (heterogenous effects in different studies) - solution: random effects 
         pooling   

      Convention:  First perform a test of homogeneity of means. 
                   If no significant inter-study variation is found, then a fixed-effects 
                   approach is adopted; otherwise a random-effects approach is adopted or study 
                   characteristics have to be identified that stratify the studies into subsets 
                   with homogeneous effects.    

      If homogeneity (null-hypothesis) is rejected, one may conclude that the study means arose 
      from two or more distinct populations. If homogeneity cannot be rejected, one would 
      conclude that the k studies share a common mean T and estimate T. The test has low power, 
      which means the probability of accepting the true alternative hypothesis (heterogeneity) is low. 
      Not rejecting the null-hypothesis is equivalent to asserting that the amount of between-study 
      variation is small.   

      If tau**2 would be known, then T could be estimated by maximum likelihood estimation (MLE). 
      In the more realistic case of unknown tau**2 it can be estimated my restricted maximum 
      likelihood (REML) or by method of moments (MOM)    */



/* Pooling macro */

%macro pooling (datein);

/* Dataset preparation */

data new; 
    set &datein;
    keep study beta var int;
	int=1;         **** Necessary as a merging parameter later on;
run;



data var(keep = value row col);    **** Needed for random effects approach later on;
   set new;
   row=_n_;        **** _n_ is SASs built-in case counter;
   col=_n_;
   value=var;      **** Specify G-matrix (variance matrix of random effects) so that (row,col) element (=diagonal) is equal to value, all other elements equal to 0;
   call symput("obs",_n_);   
run;                         /* Symput can give a value from a SAS data step to a macro variable. 
                                The value of argument2 will be given to argument1 (=name of macro variable, in this case a string) */


/* Fixed effect estimator */
/* Calculation of fixed effect (FE) estimator (maximum likelihood estimator: MLE) and test of homogeneity */

data new2;
   set new;
   w=1/var;         **** Weight: huge variance in study-specific estimate - small weight in the pooling procedure;
   beta_w=w*beta;   **** Weighted estimate;
   w2=w*w; 					
run;

proc means noprint data = new2;    **** Create sums of beta_w (s_b_w) w (s_w) and w2 (s_w2);
   var beta_w w w2;
   output out = q_data sum = s_b_w s_w s_w2;   
run;                               /* sum is an option in the means procedure
                                      and the result will be named s_b_w, s_w and s_w2 */

data q_data2;
   set q_data;
   fe=s_b_w/s_w;			  **** Fixed effects estimator after Normand;
   fe_se=sqrt(1/s_w);		  **** With standard error;
   int=1;                     **** Necessary as a merging parameter later on;
   t=fe/fe_se; 				
   prob_fe=1-probchi(t*t,1);  **** Chi-square distribution for calculating p-value for pooled effect estimate;
run;

data new3;                    **** Calculate p-value for test of homogeneity after Normand;
   merge new2 q_data2(keep = int fe);
   by int;
   wi_fe=w*(beta-fe)*(beta-fe);
run;

proc means noprint data = new3;
   var wi_fe;
   output out = q sum = q n=k;
run;                          /* sum and n are options in the means procedure
                                 and the results will be named q and k, also the output dataset will be named q */

data q2;
   set q;
   prob_hom = 1-probchi(q, k-1); **** p-value for test of homogeneity after Normand;
run;                             /* A statistical test for the homogeneity of study means is equivalent to testing:              
                                      H0: T=T1=T2=....=Tk
                                      H1: at least one Ti is different
                                    Reject H0 if p<0.05 or p<0.10 (more conservative as test has low power) */

data fixed_eff;
   merge q2(keep = prob_hom k q) q_data2(keep = fe fe_se prob_fe s_w s_b_w s_w2);  **** No by-statement necessary as only one observation per dataset;
   tau_2=max(0, (q-(k-1))/(s_w-s_w2/s_w));  **** Method of moments (MOM)-estimator for tau**2 provided by homogeneity test after DerSimonian/Laird;
   call symput ("tau_2", tau_2);
   int=1;
run;


/*  Random effect estimators */

/* Methods to be used for estimating pooled effects if heterogeneity is present:
            Method of moments (MOM)               Random Effects Estimator after DerSimonian/Laird
            Restricted maximum likelihood (REML)  Random Effects Estimator after Normand
            Maximum likelihood (ML)               Random Effects Estimator after van Houwelingen */


/* Calculation of DerSimonian/Laird estimator (MOM estimator) for the random effects model */
/* This is the estimator we are using for our exercise */
/* The DerSimonian/Laird estimator is different from the fixed effect estimator if tau_2 not equal to 0 */

data new4;
   merge fixed_eff(keep = tau_2 int) new3;
   by int;
   w_re=1/(var+tau_2);   **** Weight;
   betaw_re=beta*w_re;   **** Weighted estimate;
run;

proc means noprint data = new4;
   var betaw_re w_re;
   output out= re_dl sum = sbetaw_re sw_re n=k;;
run;

data re_dl2(keep= re_dl se_re_dl prob_re_dl);
   set re_dl;
   re_dl=sbetaw_re/sw_re;  **** MOM estimator;
   se_re_dl=sqrt(1/sw_re); **** With standard error;
   t=re_dl/se_re_dl; 				
   prob_re_dl=1-probchi(t*t,1);  **** Chi-square distribution for calculating p-value for pooled effect estimate;
run;


/* Calculation of REML estimator after Normand for the random effects model */

title "REML random effects estimator after Normand";
proc mixed data = new order = data method = reml;   **** Use order of appearance in input dataset: should be sorted by variable study;
   class study;    **** Specifies study as classification variable;
   model beta = int / s noint ddf = 10000; **** An intercept only model, print solution for fixed effects, denominator degrees of freedom high (even if small n of observations): trick for obtaining p of normal distribution;    
   random study / gdata=var s;    **** The study is specified as a random effect, print solutions for random effects, G-matrix is diagonal with study variances in diagonal;
   repeated / type=toep(1);     **** Sampling variance matrix R (variance matrix of fixed effects) is diagonal: Töplitz structure (1) causes this;
   ods output solutionF=summ;   **** Output delivery system in SAS puts solutions into a SAS datset;
run;                            **** Trick is to exchange G and R in this code as SAS does not permit fixing G, only R;

data summ2(keep= re_reml se_re_reml prob_re_reml);
   set summ;
   rename estimate = re_reml stderr = se_re_reml probt = prob_re_reml;
run;


/* Calculation of ML estimator after van Houwelingen for the random effects model */
/* This estimator converges better than the REML estimator and is less biased than the MOM estimator 
   if study-specific estimates are heterogeneous. If study-specific estimates are not heterogeneous then  
   the van Houwelingen estimator is similar to the fixed effects estimator */

data var2;
   set var;
   %do i = 1 %to &obs;   
      if row= &i then do;
      call symput("v&i",value);
      end;
   %end;
run;
%let obs1=%eval(&obs+1);    **** Evaluation of numerical statement: adding obs and 1 and writing it into obs1; 

title "ML random effects estimator after van Houwelingen";
proc mixed data = new order =data method = ml;
   class study;
   model beta =  / s cl;    **** An intercept only model (int automatically included), asking for confidence limits;
   random int / subject = study s;   **** The intercept is specified as a random effect with subject study = random study as in REML;
   repeated / group = study;    **** Each study has its own within study variance (not necessary in our case as only one estimate per study);
   parms (&tau_2) %do o=1 %to &obs; %str((&&v&o.)) %end; / eqcons=2 to &obs1;   **** Starting values of covariance parameters: the between study random component and the within-study error variances;
   ods output solutionF=summ_ml;
run;

data summ2_ml(keep= re_ml se_re_ml prob_re_ml);
   set summ_ml;
   rename estimate = re_ml stderr = se_re_ml probt = prob_re_ml;
run;


/* Create output dataset that includes all calculated estimates */

data erg;
   merge summ2 fixed_eff(keep =fe fe_se prob_fe prob_hom tau_2) re_dl2(keep=re_dl se_re_dl prob_re_dl) summ2_ml;
run;

proc print data=work.erg;
run;

data erg;
  set erg;
  daten="&datein";
run;

data all;
  retain daten;
  set erg;
  format fe fe_se re_dl prob_hom prob_fe prob_re_ml prob_re_reml re_ml re_reml se_re_dl se_re_ml se_re_reml 12.8;
run;

proc datasets library=work nolist; 
   delete erg fixed_eff new new2 new3 new4 q q2 q_data q_data2 re_dl re_dl2 summ summ2 summ2_ml summ_ml var var2; 
run; 

quit;

%mend pooling;




