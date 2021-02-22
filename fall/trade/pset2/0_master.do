*-----------------------------------------------------------*
*	 International Macro and Trade Assignment 2,
*						Tanya Rajan
* Description: This .do file runs all the Stata files for
* Questions 1, 2, and part of Question 3 on this assignment.
* 
* Note: The R and Julia code must be run separately 
*-----------------------------------------------------------*

// Setup
clear all
set more off
set maxvar 32000
timer clear
set matsize 6000
*set max_memory 16g


// Filepath 
// *** you will have to separately change filepath in the R and Julia code files
global filepath "/Users/tanyarajan/Documents/git/psets/trade/pset2"
cd ${filepath}

// Select files to run ("yes" or "no")
global runQ1 = "yes"
global runQ2 = "yes"
global runQ3 = "yes"


// Running files
do 1_pset2_clean

if "$runQ1" == "yes"{
	di "**** RUNNING QUESTION 1 CODE ****"	
	qui do 2_pset2.1.do
}

if "$runQ2" == "yes"{
	di "**** RUNNING QUESTION 2 CODE ****"	
	qui do 3_pset2.2.do
}

if "$runQ3" == "yes"{
	di "**** RUNNING QUESTION 3 CODE ****"	
	qui do 4_pset2.3.do
}

