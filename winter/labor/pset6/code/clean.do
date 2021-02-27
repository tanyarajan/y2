*--------------------------------------------
*       Cleaning NLSY Income Data
*--------------------------------------------
clear all 
global projects : env projects
global filepath ${projects}/y2/winter/labor/pset6


* Loading and naming data
insheet using ${filepath}/data/nlsy_income.csv
rename *, upper
do ${filepath}/code/nlsy_income-value-labels.do
rename *, lower

* Renaming variables of interest
rename sample_sex_1979 sex
rename caseid_1979 caseid

* Reshaping
reshape long q13_5_trunc_revised_ esr_col_, i(caseid sex) j(year)
keep q13_5_trunc_revised_* esr_col_* caseid sex year

* Renaming again
rename q13 income
rename esr employment

* Restricting sample
drop if sex == 2 // dropping females
drop if income == . | income <= 0 // missing income
drop if employment == 3 // out of labor force


* Keeping if they have at least 9 consecutive obs
xtset caseid year
gen nrun = .
by caseid: replace nrun = cond(L.nrun == ., 1, L.nrun + 1)
by caseid: egen maxrun = max(nrun)
keep if maxrun > 9 

* Removing obs
gen runcnt = L.nrun
drop if runcnt == .

* Lags
gen lag1 = L.income
gen fd_inc = income - L.income
drop if fd_inc == .

* Save
outsheet using ${filepath}/data/nlsy_clean.csv, comma nolabel replace
