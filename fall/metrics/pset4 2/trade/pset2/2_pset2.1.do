*-----------------------------------------------------------*
*	 International Macro and Trade Assignment 2, Table 1
*						Tanya Rajan
* Description: This file uses data associated with Head
* and Mayer (2014) to estimate gravity equations using 
* various Stata fixed effect estimation strategies. The
* analysis uses only years 2000-2006.
*-----------------------------------------------------------*
// Setup
clear all
set more off
set maxvar 32000
timer clear
set matsize 6000
*set max_memory 16g


// filepath set in 0_master.do
use ${filepath}/data/clean.dta, clear


************************************
*** 	Table 1 Regressions		 ***
************************************
// Restrictions
keep if year >= 2000 & year <= 2006 
global subset1 "if flow != 0"
global subset2 "if flow != 0 & year == 2006" 

// Singleton observations when these variables = 1
egen single_exp = count(flow), by(expyear)
egen single_imp = count(flow), by(impyear)

*** Using reg ***
timer clear 1
timer on 1
qui eststo m1: reg lnflow i.exporter#i.year i.importer#i.year ///
			lndist contig comlang_off ${subset1}
timer off 1 
timer list 1
estadd local time = r(t1), replace


// Smaller subset (.5210 secs)
timer clear 20
timer on 20
qui reg lnflow i.exporter#i.year i.importer#i.year ///
			lndist contig comlang_off ${subset2}
timer off 20
timer list 20




*** Using xtreg ***
xtset expyear
timer clear 2
timer on 2
qui eststo m2: xtreg lnflow i.impyear lndist contig comlang_off ${subset1}, fe 
timer off 2
timer list 2
estadd local time = r(t2), replace

// Smaller subset (1.349 secs)
timer clear 21
timer on 21
qui xtreg lnflow i.impyear lndist contig comlang_off ${subset2}, fe 
timer off 21
timer list 21



*** Using areg ***
timer clear 4
timer on 3
qui eststo m3: areg lnflow i.impyear lndist contig comlang_off ///
			${subset1}, absorb(expyear)
timer off 3
timer list 3
estadd local time = r(t3), replace

// Smaller subset (0.279 secs)
timer clear 22
timer on 22
qui areg lnflow i.impyear lndist contig comlang_off ///
			${subset2}, absorb(expyear)
timer off 22
timer list 22



*** Using hdfe ***
timer clear 4
timer on 4
qui eststo m4: reghdfe lnflow lndist contig comlang_off ${subset1}, ///
			absorb(expyear impyear) 
timer off 4
timer list 4
estadd local time = r(t4), replace

// Smaller subset (0.1640)
timer clear 23
timer on 23
qui reghdfe lnflow lndist contig comlang_off ${subset2}, ///
			absorb(expyear impyear) 
timer off 23
timer list 23
			


*** Writing out to Latex ***
esttab m1 m2 m3 m4 using ${filepath}/tables/tr_table1.tex, se replace label ///
		keep(lndist) s(time N r2 r2_a, label("Time" "N" "$ R^2$" "Adj. $ R^2$")) ///
		nostar compress ///
		mtitles("\textbf{reg}" "\textbf{xtreg}" "\textbf{areg}" "\textbf{reghdfe}") 

			
			
				
				
