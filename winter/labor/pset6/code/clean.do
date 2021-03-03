*--------------------------------------------
*       Cleaning Data to Replicate the
* 				Hryshko Analysis
*--------------------------------------------
* Note: ssc install tsspell if you haven't already
* also may need to install carry forward 
ssc install tsspell
ssc install carryforward

clear all 
global projects : env projects
global filepath ${projects}/psets_heckmanlabor/pset6
*--------------------------------------------
* 				Cleaning CPI data
*-------------------------------------------- 
insheet using ${filepath}/data/cpi.csv

* Averaging per year
collapse (mean) value, by(year)
rename value cpi
sort year

* Normalizing 1982 to 100 as in the paper
local base = cpi
replace cpi = cpi/`base'

save ${filepath}/data/cpi_adjustment.dta, replace


*--------------------------------------------
* 				Cleaning NLSY data
*-------------------------------------------- 
* Loading and naming data
insheet using ${filepath}/data/nlsy_income.csv, clear
rename *, upper
do ${filepath}/code/nlsy_income-value-labels.do
rename *, lower

* Renaming variables of interest
rename sample_sex_1979 sex
rename caseid_1979 caseid
rename fam_1b_1979 age

* Reshaping
reshape long q13_5_trunc_revised_ esr_col_ q3_4_ q3_3_, i(caseid sex) j(year)
keep q13_5_trunc_revised_* esr_col_* q3_4_* q3_3* caseid sex year sampweight age

* Renaming again
rename q13 income
rename esr employment
rename q3_4 edu_complete
rename q3_3 edu_attend
rename sampweight weights

***** Age and Education Variables *****

* Restricting sample
drop if sex == 2 // dropping females
replace income = . if income <= 0 // missing income
replace income = . if employment == 3 // out of labor force

* Generating college dummy
replace edu_complete = . if edu_complete < 0
replace edu_attend = . if edu_attend < 0
sort caseid year
bys caseid: carryforward edu_complete, gen(education) // imputing education for future years
egen nonmiss_edu = count(edu_complete), by(caseid)
drop if nonmiss_edu == 0

// College dummy
gen college = (education > 12)
replace college = . if education == .

* Generating age polynomials
replace age = (year - 1979) + age
drop if age <= 25 | age > 64
label val age
gen age2 = age^2
gen age3 = age^3

***** Keeping Every Other Period *****
* generating period indicator to deal with 2-year gaps that start in 1996 
* note the start year is 1982 bc the standardized income variable starts then
drop if inlist(year, 1979, 1981, 1983, 1985, 1987, 1989, 1991, 1993, 1995)

gen period = year - 1982
global yrlist 1984 1986 1988 1990 1992 1994 1996 1998 2000 ///
			2002 2004 2006 2008 2010 2012 2014 2016
global num : list sizeof global(yrlist)
forval y = 1/$num{
	local yr : word `y' of ${yrlist}
	replace period = `y' if year == `yr'
}
xtset caseid period


***** Income Variables *****

* Adjusting by CPI
rename income unadj_inc
merge m:1 year using ${filepath}/data/cpi_adjustment.dta
drop if _merge == 2
drop _merge
sort caseid period
gen income = unadj_inc / cpi

* Restricting if real change in income too high or low as in paper
gen inc_change = (income - L.income)/income * 100
replace income = . if inc_change > 500 | inc_change < -80 


* Log income
sort caseid period
gen lninc = ln(income)
gen d_lninc = lninc - L.lninc
egen ln_var = sd(lninc), by(year)
replace ln_var = ln_var^2
preserve
collapse (mean) ln_var, by(year)
*twoway line ln_var year if year >= 1984, graphregion(color(white))
restore

* Residual income: regression of difference in log earnings on edu/age 
reg d_lninc age age2 age3 college age#college age2#college age3#college
predict y, residuals

* Marking consecutive runs of data
drop if year < 1984
gen nrun = .
by caseid: replace nrun = cond(L.nrun == . | y == . | L.y== ., 1, L.nrun + 1)
by caseid: egen maxrun = max(nrun)

* Keeping only longest consecutive spell (and only if 9 consecutive observations)
tsspell income, f(nrun==1)
egen spellmax = count(_spell), by(caseid _spell)
keep if spellmax == maxrun
keep if spellmax > 9 

*--------------------------------------------
* 		Setting up lags for GMM
*-------------------------------------------- 
keep caseid period maxrun y
sum maxrun
global maxT = r(max) 

* T 
sort caseid period
gen T = .
bys caseid: replace T = _n

drop period maxrun 
order caseid T y

* Lags needed
xtset caseid T
gen ly1 = L.y
gen ly2 = L2.y

global loopT = $maxT - 1
forval t=1/$loopT{
	gen fy`t' = F`t'.y
}

* Save for analysis in R
outsheet using ${filepath}/data/nlsy_clean.csv, comma nolabel replace
