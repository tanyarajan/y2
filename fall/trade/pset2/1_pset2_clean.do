*-----------------------------------------------------------*
*	 International Macro and Trade Assignment 2,
*						Tanya Rajan
* Description: This .do file cleans the data from Head
* and Mayer (2014) and creates necessary variables to run
* the gravity estimation for this assignment.
*-----------------------------------------------------------*

// Setup
clear all
set more off
set maxvar 32000
timer clear
set matsize 6000
*set max_memory 16g


************************************
*** 	Cleaning Variables 		 ***
************************************
// Load data (filepath set in 0_master.do)
use flow distw iso_o iso_d year contig comlang_off gdp_o gdp_d ///
		using ${filepath}/data/col_regfile09.dta, clear

// Log variables
gen lnflow = ln(flow)
gen lnflow1 = ln(flow + 1)
gen lndist = ln(distw) 
gen lnexpGDP = ln(gdp_o)
gen lnimpGDP = ln(gdp_d)

// Generating dummies
encode(iso_o), gen(exporter)
encode(iso_d), gen(importer)

egen bilat_id = group(iso_o iso_d)
egen expyear = group(exporter year)
egen impyear = group(importer year)


// Labels
label var lndist "Log Distance"
label var lnflow1 "Log(Flow + 1)"
label var lnflow "Log Flow"
label var contig "Contiguity"
label var comlang_off "Common Language"

// Reducing size of data
drop iso_o iso_d 
compress _all

// Saving cleaned data for R and Julia estimations
save ${filepath}/data/clean.dta, replace


