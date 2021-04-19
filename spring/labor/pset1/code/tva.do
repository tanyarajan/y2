*--------------------------------------------
*       Teacher Value Added Estimates
*--------------------------------------------
* Note: ssc install the following: 
ssc install unique

* Setup
clear all 
global projects : env projects
global filepath ${projects}/Education/pset1
cd ${filepath}

use data/va-data.dta, clear

*--------------------------------------------
* 	  Constructing the Control Vector
*-------------------------------------------- 
* dataset is identified at the student year level
xtset id_student year

// Renaming
*ren lag_score_std lagscore
ren score_std scores
gen lagscore = L.scores
ren id_grade grade

// Lag individual test scores
gen lagscore2 = lagscore^2
gen lagscore3 = lagscore^3
gen lagscore_grade = lagscore*grade
gen lagscore2_grade = lagscore2*grade
gen lagscore3_grade = lagscore3*grade

// Identifying classes
* each class is identified at id_class, id_school, id_grade level
* class means this is a 4th grade class with ID #1 at school #X
egen class = group(id_class id_school year grade)
egen csize = count(score), by(class year)

// Class test scores
sort id_student year
egen class_score = mean(score), by(class year)
gen cs_lag = L.class_score
gen cs_lag2 = cs_lag^2
gen cs_lag3 = cs_lag^3

// School test scores
egen school_score = mean(score), by(id_school year)
gen ss_lag = L.school_score
gen ss_lag2 = ss_lag^2
gen ss_lag3 = ss_lag^3

// Grade and year dummies
tab grade, gen(g_)
tab year, gen(y_)

// Class and school-year avgs of demographic variables
egen cavg_inc = mean(hh_income), by(class)
egen cavg_medu = mean(m_education), by(class)
egen syavg_inc = mean(hh_income), by(id_school year)
egen syavg_medu = mean(m_education), by(id_school year)


*--------------------------------------------
* 	  Step 1: Residualizing Test Scores
*-------------------------------------------- 

// Regression on X vector and teacher fixed effects
global xvec lagscore lagscore2 lagscore3 c.lagscore#i.grade c.lagscore2#i.grade ///
			c.lagscore3#i.grade cs_lag cs_lag2 cs_lag3 c.cs_lag#i.grade c.cs_lag2#i.grade ///
			c.cs_lag3#i.grade ss_lag ss_lag2 ss_lag3 c.ss_lag#i.grade c.ss_lag2#i.grade ///
			c.ss_lag3#i.grade cavg_inc cavg_medu syavg_inc syavg_medu csize ///
			i.grade i.year
			
reghdfe scores ${xvec}, absorb(teach_eff = i.id_teacher) resid(residual) nocons
global mse = e(rmse)^2


// Predicting (leaving out teacher fixed effects)
gen A = residual + teach_eff


*--------------------------------------------
* 	Step 2: Autocovariance of Test Scores
*--------------------------------------------
// Computing the degrees of freedom correction outlined in Appendix
unique id_student
global Nval = r(unique)
unique class
global Cval = r(unique)
global Kval : list sizeof global(xvec)
global df_correction = ((${Nval} - 1)/(${Nval} - ${Kval} - ${Cval} + 1))

// Calculating sigma_epsilon (not sure whether to use mse from regression or 
// the one computed below)
gen residual2 = A^2
egen ssr = total(residual2)
replace ssr = ssr*(1/${Nval})
global s_ep = ssr*${df_correction}
		
// Calculating sigma_A
corr A, cov
global s_A = r(Var_1)*${df_correction}


/* Can't calculate sigma_A0 because we don't have data with multiple
 classes per teacher in a year. But also maybe don't need this whole section
 because we don't need to precision weight things. */
 
 
// Creating class means of A
egen class_A = mean(A), by(class)


// Saving data
save data/va_lean.dta, replace
outsheet using data/va_clean.csv, comma nolab replace



 
// Collapsing to the teacher-year level
collapse (mean) A (count) id_student, by(id_teacher year class grade)
ren id_student ns
xtset id_teacher year
drop if A == .


// Generating lags for autocov
sort id_teacher year
forval i = 1/6{
	gen A_l`i' = L`i'.A
	gen ns_l`i' = L`i'.ns
	gen l`i'_weight = ns + ns_l`i'
}

// Calculating autocovs (weighted by sum of students in pair of classrooms)
forval i = 1/6{
	qui corr A A_l`i' [aw=l`i'_weight], cov
	global s_A`i' = r(cov_12)
	qui corr A A_l`i' [aw=l`i'_weight]
	global r_`i' = r(rho)
	di "Lag `i', autocov: " ${s_A`i'} " and autocorr: " ${r_`i'}
}

/* Something is wrong here: using their vector of controls X as defined on pg 2605,
 I get autocorrelations that are an order of magnitude smaller than what's 
 reported in fig1/tab2 of the paper. However if i include hh_income in the x vector,
 I get autocorrelations much closer to what they report... not sure what to make of this.
 We are also missing a lot of the data they use to control in that X vector. */


*--------------------------------------------
* 		Step 3: Estimating Teacher VA
*--------------------------------------------

/* might be better to do this in R or other programming language that is more
matrix friendly */

save data/va_collpaseclean, replace
outsheet using data/va_collapseclean.csv, comma nolab replace
			
