*-----------------------------------------------------------*
*	 International Macro and Trade Assignment 2, Table 2
*						Tanya Rajan
* Description: This file uses data associated with Head
* and Mayer (2014) to estimate various specifications of
* the gravity equation including log-linear and CES versions
* both including and excluding observations with 0 flow. The
* analysis uses only years 2004-2006.
*-----------------------------------------------------------*
// Setup
clear all
set more off
set maxvar 32000
timer clear
set matsize 6000


// filepath set in 0_master.do
use ${filepath}/data/clean.dta, clear


************************************
*** 	Table 2 Regressions		 ***
************************************
// restricting set
drop if year < 2004

// subset with no 0s
global subset1 "if flow != 0"

*** Log linear ***
timer clear 1
timer on 1
qui eststo s1: reghdfe lnflow lndist contig comlang_off ${subset1}, ///
			absorb(expyear impyear) resid
predict resid, resid
timer off 1
timer list 1
estadd local time = r(t1), replace
estadd local zeros "No"
estadd local flow1 "No"
estadd local regtype "Log-Linear"


*** Log linear with log(flow+1) ***
timer clear 2
timer on 2
qui eststo s2: reghdfe lnflow1 lndist contig comlang_off ${subset1}, ///
			absorb(expyear impyear)
timer off 2
timer list 2
estadd local time = r(t2), replace
estadd local zeros "No"
estadd local flow1 "Yes"
estadd local regtype "Log-Linear"


*** Log linear with log(flow+1), include 0 ***
timer clear 3
timer on 3
qui eststo s3: reghdfe lnflow1 lndist contig comlang_off, ///
			absorb(expyear impyear)
timer off 3
timer list 3
estadd local time = round(r(t3), .001), replace
estadd local zeros "Yes"
estadd local flow1 "Yes"
estadd local regtype "Log-Linear"

				
*** CES specification with ppml ***
// Generating dummies since ppml doesn't allow factor variables
tab expyear, gen(ey_)
tab impyear, gen(iy_)

timer clear 4
timer on 4
qui eststo s4: ppml flow lndist ey_* iy_* comlang_off contig
timer off 4
timer list 4
estadd local time = r(t4), replace
estadd local zeros "Yes"
estadd local flow1 "No"
estadd local regtype "Poisson"


*** CES specification with poi2hdfe ***
timer clear 5
timer on 5
qui eststo s5: poi2hdfe flow lndist contig comlang_off, id1(expyear) id2(impyear)
timer off 5
timer list 5
estadd local time = r(t5), replace
estadd local zeros "Yes"
estadd local flow1 "No"
estadd local regtype "Poisson"


				
*** CES specification with ppml_panel_sg ***
timer clear 6
timer on 6
qui eststo s6: ppml_panel_sg flow lndist contig comlang_off, ///
	exporter(exporter) importer(importer) year(year) nopair olsguess nocheck
timer off 6
timer list 6
estadd local time = r(t6), replace
estadd local zeros "Yes"
estadd local flow1 "No"
estadd local regtype "Poisson"



*** CES specification with ppmlhdfe ***
timer clear 7
timer on 7
qui eststo s7: ppmlhdfe flow lndist contig comlang_off, ///
	absorb(expyear impyear)
timer off 7
timer list 7
estadd local time = r(t7), replace			
estadd local zeros "Yes"	
estadd local flow1 "No"		
estadd local regtype "Poisson"


*** CES specification with ppmlhdfe exclude 0 ***
timer clear 8
timer on 8
qui eststo s8: ppmlhdfe flow lndist contig comlang_off ${subset1}, ///
	absorb(expyear impyear)
timer off 8
timer list 8
estadd local time = r(t8), replace		
estadd local zeros "No"
estadd local flow1 "No"
estadd local regtype "Poisson"


************************************
***   Heteroskedasticity Tests	 ***
************************************

*** Breusch Pagan Test ***
// without zeroes
qui reg lnflow i.exporter#i.year i.importer#i.year ///
			lndist contig comlang_off ${subset1}
hettest

// with zeroes
qui reg lnflow i.exporter#i.year i.importer#i.year ///
			lndist contig comlang_off
hettest
							
				
				
*** Writing out to Latex ***
esttab s1 s2 s3 s4 s5 s6 s7 s8 using ${filepath}/tables/tr_table2.tex, ///
		se replace label keep(lndist contig comlang_off) ///
		s(N r2 time zeros regtype flow1, ///
		label( "N" "$ R^2$" "Time" "Include Zeroes?" "Method" "log(Flow + 1)?") ///
		fmt(%10.0fc %10.3fc %10.3fc)) ///
		nostar compress ///
		mtitles("\textbf{reghdfe}" "\textbf{reghdfe}" "\textbf{reghdfe}" ///
		"\textbf{ppml}" "\textbf{poi2hdfe}" "\textbf{ppml_panel_sg}" ///
		"\textbf{ppmlhdfe}" "\textbf{ppmlhdfe}")
			

*** Graph residuals ***
scatter resid lndist, graphregion(color(white)) msize(.3) mcolor(dkorange)
graph export ${filepath}/figures/tr_residuals.png, replace

