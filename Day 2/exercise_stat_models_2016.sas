**********************************************************************
Project:      Advanced Epi / Env & Occ Epi WS15/16
Program:      exercises_stat_models.sas
Created:      18.01.2016
Data sets:    mort, mort_cc, escape								
**********************************************************************;

* create library "class" *;

*LIBNAME class "D:\0_toCopy\Day04_Statistical_Methods\Exercise\data";
LIBNAME class "/home/msc1426/Dokumente/EnvEpi-Scripts/Day 2";

/***************************************/
/*1. CVD mortality - POISSON REGRESSION*/
/***************************************/

* Read in the dataset "mort" in the temporary work space *;
DATA mort;
SET class.mort;
	log_pop=log(pop_augs); /*logarithm of Augsburg population*/
	trend=trend/1000; /*divide by 1000 for larger effect and std.error estimates*/
	temp2=M24T_A_L04**2; /*quadratic and cubic term of temperature*/
	temp3=M24T_A_L04**3;
RUN;


* Description of number of cases per day *;
PROC MEANS DATA=mort MAXDEC=2 MEAN VAR MIN MAX;
	VAR totmort_augst cardiomort_augst cardiomort_augst_sex0 cardiomort_augst_sex1;
RUN;


* Plot of daily cardiovascular mortality during the study period *; 
PROC GPLOT DATA=mort;
title "Plot 1 sm=20";
SYMBOL1 V=POINT C=black;/*symbol for daily dots*/
SYMBOL2 V=NONE I=sm20 C=red; /*spline interpolation ->smooth curve*/
	PLOT cardiomort_augst*date=1 cardiomort_augst*date=2/OVERLAY;
RUN;QUIT;

PROC GPLOT DATA=mort;
title "Plot 1 sm=05";
SYMBOL1 V=POINT C=black;/*symbol for daily dots*/
SYMBOL2 V=NONE I=sm05 C=blue; /*spline interpolation ->smooth curve*/
	PLOT cardiomort_augst*date=1 cardiomort_augst*date=2/OVERLAY;
RUN;QUIT;

PROC GPLOT DATA=mort;
title "Plot 1 sm=99";
SYMBOL1 V=POINT C=black;/*symbol for daily dots*/
SYMBOL2 V=NONE I=sm99 C=green; /*spline interpolation ->smooth curve*/
	PLOT cardiomort_augst*date=1 cardiomort_augst*date=2/OVERLAY;
RUN;QUIT;


* Seasonal variation of mortality in 2000 *;
PROC GPLOT DATA=mort;
title "Seasonal variation of mortality in year 2000";
SYMBOL1 V=dot C=black;
SYMBOL2 V=NONE I=sm50 C=red; 
AXIS1 ORDER=('01jan00'd to '01jan01'd by month);
	PLOT cardiomort_augst*date=1 cardiomort_augst*date=2/OVERLAY HAXIS=axis1;
RUN;QUIT;

PROC GPLOT DATA=mort;
title "Seasonal variation of mortality in year 2001";
SYMBOL1 V=dot C=black;
SYMBOL2 V=NONE I=sm50 C=red; 
AXIS1 ORDER=('01jan01'd to '01jan02'd by month);
	PLOT cardiomort_augst*date=1 cardiomort_augst*date=2/OVERLAY HAXIS=axis1;
RUN;QUIT;

PROC GPLOT DATA=mort;
title "Seasonal variation of mortality in year 1990";
SYMBOL1 V=dot C=black;
SYMBOL2 V=NONE I=sm50 C=red; 
AXIS1 ORDER=('01jan90'd to '01jan91'd by month);
	PLOT cardiomort_augst*date=1 cardiomort_augst*date=2/OVERLAY HAXIS=axis1;
RUN;QUIT;


* Do holidays (HOLIDAY, 1=yes) and weekdays (DOW, 1=Sunday) influence the number of daily mortality cases? Use PROC MEANS *;
* Cardiovascular mortality by holiday *;
PROC MEANS DATA=mort MAXDEC=2 MEAN STD MIN MAX ;
title "Cardiovascular mortality by holiday";
	VAR cardiomort_augst;
	CLASS holiday;
RUN;

* Cardiovascular mortality by day of the week *;
PROC MEANS DATA=mort MAXDEC=2 MEAN STD MIN MAX ;
title "Cardiovascular mortality by DOW";
	VAR cardiomort_augst;
	CLASS dow;
RUN;


* Poisson regression: Effect of PM10 on cardiovascular mortality *;
PROC GENMOD DATA = mort;
title "Effect of PM10 on cardiovascular mortality";
	CLASS  season dow (REF=first)/PARAM=ref;
	MODEL cardiomort_augst=M24PM10_A_L0 M24T_A_L04 temp2 temp3 M24RH_A_L04 M24BP_A_L04 trend season dow / DIST=poisson LINK=log;
RUN;

* take overdispersion into account *;
PROC GENMOD DATA = mort;
title "Effect of PM10 on cardiovascular mortality w/overdispersion";
	CLASS  season dow(REF=first)/PARAM=ref;
	MODEL cardiomort_augst=M24PM10_A_L0 M24T_A_L04 temp2 temp3 M24RH_A_L04 M24BP_A_L04 trend season dow/ DIST=poisson LINK=log SCALE=pearson ;
RUN;

* take changing population into account *;
PROC GENMOD DATA = mort;
title "Effect of PM10 on cardiovascular mortality w/changing population";
	CLASS  season dow(REF=first)/PARAM=ref;
	MODEL cardiomort_augst=M24PM10_A_L0 M24T_A_L04 temp2 temp3 M24RH_A_L04 M24BP_A_L04 trend season dow/ DIST=poisson LINK=log SCALE=pearson OFFSET=log_pop;
	ESTIMATE 'PM10' M24PM10_A_L0 10/EXP; /*calculation of RR, label ('PM10') must be specified
										increment must be specified, here=10*/
RUN;


* Which lag of NO2 has the strongest effect on cardiovascular mortality? *;
%MACRO lags (poll=);
	%DO j=0 %TO 5;
	%LET i=&j;
	%IF &j=5 %THEN %DO;
	%LET i=04;
	%END;
		PROC GENMOD DATA = mort;
			title "Effect of NO2 on CV mortality";
			CLASS  season dow(REF=first)/PARAM=ref;
			MODEL cardiomort_augst=M24&poll._A_L&i M24T_A_L04 temp2 temp3 M24RH_A_L04 M24BP_A_L04 trend season dow/ DIST=poisson LINK=log SCALE=pearson OFFSET=log_pop;
			ESTIMATE "&poll." M24&poll._A_L&i 10/EXP; 
		RUN;
	%END;
%MEND;

%lags (poll=NO2);

* Which lag of NO2 has the strongest effect on cardiovascular mortality? *;
%MACRO lags (poll=);
	%DO j=0 %TO 5;
	%LET i=&j;
	%IF &j=5 %THEN %DO;
	%LET i=04;
	%END;
		PROC GENMOD DATA = mort;
			title "Effect of NO2 on CV mortality";
			CLASS  season dow(REF=first)/PARAM=ref;
			MODEL cardiomort_augst=M24&poll._A_L&i M24T_A_L04 temp2 temp3 M24RH_A_L04 M24BP_A_L04 trend season dow/ DIST=poisson LINK=log SCALE=pearson OFFSET=log_pop;
			ESTIMATE "&poll." M24&poll._A_L&i 10/EXP; 
		RUN;
	%END;
%MEND;

%lags (poll=NO2);

* Which lag of NO2 has the strongest effect on cardiovascular mortality? *;
%MACRO lags (poll=);
	%DO j=0 %TO 5;
	%LET i=&j;
	%IF &j=5 %THEN %DO;
	%LET i=04;
	%END;
		PROC GENMOD DATA = mort;
			title "Effect of NO2 on CV mortality";
			CLASS  season dow(REF=first)/PARAM=ref;
			MODEL cardiomort_augst=M24&poll._A_L&i M24T_A_L04 temp2 temp3 M24RH_A_L04 M24BP_A_L04 trend season dow/ DIST=poisson LINK=log SCALE=pearson OFFSET=log_pop;
			ESTIMATE "&poll." M24&poll._A_L&i 10/EXP; 
		RUN;
	%END;
%MEND;

%lags (poll=NO2);





/****************************************************/
/*2. CVD mortaltiy - CONDITIONAL LOGISTIC REGRESSION*/
/****************************************************/

* Read in the data set *;
DATA mort_cc;
	SET class.mort_cc;
	temp2=M24Temp_Lag04**2; /*quadratic and cubic term of temperature*/
	temp3=M24Temp_Lag04**3;
	time=2-cvd; /*artifical variable for time of death - 1 (case), 2 (control, later/other time-point)*/
RUN;


* Cond. log. reg.: Effect of PM10 on cardiovascular mortality *;
PROC LOGISTIC DATA=mort_cc;
title "Cond. log. reg.: Effect of PM10 on cardiovascular mortality";
	STRATA id; /*this statement is only available since SAS version 9*/
	MODEL cvd(EVENT='1')= M24PM10_Lag0 M24Temp_Lag04 temp2 temp3 M24RH_Lag04 M24BP_Lag04 trend fluepi;
RUN;

PROC PHREG DATA=mort_cc;
	STRATA id;
	MODEL time*cvd(0) = M24PM10_Lag0 M24Temp_Lag04 temp2 temp3 M24RH_Lag04 M24BP_Lag04 trend fluepi/TIES=discrete;
RUN;


* Which lag of NO2 has the strongest effect on cardiovascular mortality? *;
%MACRO lags2 (poll=);
	%DO j=0 %TO 5;
	%LET i=&j;
	%IF &j=5 %THEN %DO;
	%LET i=04;
	%END;
		PROC LOGISTIC DATA = mort_cc;
		title "Which lag of NO2 has the strongest effect on cardiovascular mortality?";
			STRATA id;
			MODEL cvd(EVENT='1')= M24&poll._Lag&i. M24Temp_Lag04 temp2 temp3 M24RH_Lag04 M24BP_Lag04 trend fluepi;
		RUN;
	%END;
%MEND;

%lags2 (poll=NO2);

* Explain effects of the covariates *;
PROC LOGISTIC DATA=mort_cc;
title "Complete: Cond. log. reg.: Effect of NO2 on cardiovascular mortality";
	STRATA id; /*this statement is only available since SAS version 9*/
	MODEL cvd(EVENT='1')= M24PM10_Lag0 M24Temp_Lag04 temp2 temp3 M24RH_Lag04 M24BP_Lag04 trend fluepi;
RUN;

PROC LOGISTIC DATA=mort_cc;
title "Lag04 only: Cond. log. reg.: Effect of NO2 on cardiovascular mortality";
	STRATA id; /*this statement is only available since SAS version 9*/
	MODEL cvd(EVENT='1')= M24PM10_Lag0 M24Temp_Lag04 temp2 temp3 M24RH_Lag04 M24BP_Lag04 trend fluepi;
RUN;






* Rerun the model after exclusion of barometric pressure *;
PROC LOGISTIC DATA=mort_cc;
title "Complete: Cond. log. reg.: Effect of NO2 on cardiovascular mortality";
	STRATA id; /*this statement is only available since SAS version 9*/
	MODEL cvd(EVENT='1')= M24PM10_Lag0 M24Temp_Lag04 temp2 temp3 M24RH_Lag04 M24BP_Lag04 trend fluepi;
RUN;






* How does the AIC change and what does that mean with respect to the model fit? *;





/***************************************/
/*3. CVD mortaltiy - SURVIVAL ANALYSIS*/
/***************************************/

* read data of the ESCAPE study *;
DATA escape;
SET class.escape;
RUN;


* Description of the Outcome variables *;
* How many individuals died during the follow up because of CVD (cvdmort_w6a2)? 
* How long is the average follow up time (fu_time)? *;
* Perform the description for the two surveys (study) separately *;
PROC FREQ DATA=escape;
	TABLE cvdmort_w6a2*study;
RUN;

PROC MEANS DATA=escape;
	VAR fu_time;
	CLASS study;
RUN;


* Kaplan Meier estimate by sex *;
PROC LIFETEST DATA=escape PLOTS=(s,ls, lls) CS=plus ES=none ;
	TIME fu_time*cvdmort_w6a2(0);
	STRATA sex;
	WHERE study="S3";/*only S3 participants*/
RUN;
* What is the probability to survive at least 5371 days for women? *;
*0,9617;

* Further options *;
PROC LIFETEST DATA=escape CS=none ES=none OUTSURV=surv TIMELIST=(1000,2500,5000) REDUCEOUT;
	TIME fu_time*cvdmort_w6a2(0);
RUN;

* What is the probability to survive at least 1000, 2500, and 5000 days? *;
* 1000 = 0.9952
  2500 = 0.9857
  5000 = 0.9579;

* How many persons died/dropped out at this time points? *;
* 1000 n = 41
  2500 n = 122
  5000 n = 266;

* Cox regression: effects of covariates and PM10 on mortality *;
PROC PHREG DATA=escape;
	CLASS study(ref="S3") sex (ref="0") smoking (ref=last) fruit (ref=first) 
	rawveg (ref=first) mar_stat (ref="2") edulev (ref=first) etsa 
	(ref=first) etsb  (ref=first) occstatus (ref=first) /param=ref;
	MODEL fu_time*cvdmort_w6a2(0)=age study sex smoking ncig_lif_w6 
	smoke_yrs_w6 fruit rawveg alc_current alc_current2 bmi bmi2 mar_stat
	edulev etsa etsb  occstatus gc_linc5000_07 pm10/RL;
RUN;


* Age as time variable *;
* How do effects of edulev and PM10 change? *;
PROC PHREG DATA=escape;
	CLASS study(ref="S3") sex (ref="0") smoking (ref=last) fruit (ref=first) 
	rawveg (ref=first) mar_stat (ref="2") edulev (ref=first) etsa 
	(ref=first) etsb  (ref=first) occstatus (ref=first) /param=ref;
	MODEL (age_visit,age_endfup)*cvdmort_w6a2(0)=study sex smoking ncig_lif_w6 
	smoke_yrs_w6 fruit rawveg alc_current alc_current2 bmi bmi2 mar_stat
	edulev etsa etsb  occstatus gc_linc5000_07 pm10/RL;
RUN;


* Graphical assessment of PH assumption *;
PROC LIFETEST DATA=escape PLOTS=(lls) CS=none ES=none ;
	TIME fu_time*cvdmort_w6a2(0);
	STRATA mar_stat; /* marital status */
RUN;

* Assessment of PH assumption - Test *;
PROC PHREG DATA=escape;
	CLASS mar_stat (ref="2");
	MODEL fu_time*cvdmort_w6a2(0)=mar_stat ms1_time ms3_time ms4_time;
	ms1_time=(mar_stat=1)*log(fu_time);/*create new variables*/
	ms3_time=(mar_stat=3)*log(fu_time);
	ms4_time=(mar_stat=4)*log(fu_time);
	PH_test: TEST ms1_time=0, ms3_time=0, ms4_time=0; 
RUN;

* Check the PH assumption for other covariates *;

/*NOT IN STUDENT VERSION*/
PROC LIFETEST DATA=escape PLOTS=(lls) CS=none ES=none noprint;
	TIME fu_time*cvdmort_w6a2(0);
	STRATA etsb; 
RUN;

/*NOT IN STUDENT VERSION*/
* Assessment of PH assumption - Test - smoking*;
PROC PHREG DATA=escape;
	CLASS smoking (ref=last);
	MODEL fu_time*cvdmort_w6a2(0)= smoking sm1_time sm2_time ;
	sm1_time=( smoking =1)*log(fu_time);/*create new variables*/
	sm2_time=( smoking =2)*log(fu_time);
PH_test: TEST sm1_time=0, sm2_time; 
RUN;
