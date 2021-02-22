*-----------------------------------------------------------*
*	 International Macro and Trade Assignment 2, Table 3
*						Tanya Rajan
* Description: This file uses data associated with Head
* and Mayer (2014) to run the fixed-effects log-linear version
* of the gravity regression on all years from 1948-2006.
*-----------------------------------------------------------*
// Setup
clear all
set more off
set maxvar 32000
timer clear
set matsize 6000
version 12.0

// filepath set in 0_master.do
use ${filepath}/data/clean.dta, clear

************************************
*** 	Table 1 Regressions		 ***
************************************
// Restrictions
global subset1 "if flow != 0"


*** Using reg ***
timer clear 1
timer on 1
qui eststo m1: reghdfe lnflow lndist contig comlang_off ${subset1}, ///
			absorb(expyear impyear) vce(robust)
timer off 1 
timer list 1
estadd local time = r(t1), replace

*** Writing out to Latex ***
esttab m1 using ${filepath}/tables/tr_table3_st.tex, se replace label ///
		keep(lndist contig comlang_off) s(time N, label("Time" "N")) nostar mtitle("reg") compress

			
			
				
				
