***************************************
*     Cleaning Data and Producing     
* Tables 1, 2, 5 and Figure 1 from 
*           Asker (2007)
***************************************


clear all 
global projects : env projects
global filepath ${projects}/y2/winter/io/pset3

cd $filepath

*--------------------------------------------
* 				Cleaning CPI data
*-------------------------------------------- 
import delimited using data/ps3.csv, varnames(2)

// rename things
rename datanueer id
rename realisationinfinalauation target_price
rename estdatedmindum est_min_dum
rename estdatedmaxdum est_max_dum
rename aatalogpriae cat_price
rename estdatemin est_min
rename estdatemax est_max 
rename exalusivelyus usonly


// ring wins 
egen highbid = max(bid), by(lot house date)
egen medbid = median(bid), by(lot house date)
egen mnbid = mean(bid), by(lot house date)
gen ringwin = (target_price <= highbid)

// unqiue lots
egen uniquelots = tag(lot date house)
replace uniquelots = . if uniquelots == .
egen auction = group(lot date house)


// value and profit
gen value = (est_min + est_max)/2  if uniquelots == 1
gen ringvalue = value if ringwin == 1 & uniquelots == 1
egen totprofit = total(profit), by(lot house date)

// number of bdiders and wins by bidder
egen nbidders = count(bidder), by(lot house date)
gen bidder_win = (rank == 1)
gen bidder_win2 = (rank == 1 & nbidders >= 2)
gen auction2 = auction if nbidders >=2

// side payments
gen sidepay = (netpayment > 0) if nbidders >= 2
gen payside = (netpayment < 0) if nbidders >= 2


// all manipulations
save data/allvars, replace


*** Keep only data needed for analysis ***

// Data on knockout bids
preserve
keep if nbidders == 2

tab bidder, gen(rm_)
sort auction
tab auction, gen(auc_)
keep bid bidder auction auc_* rm_* est_min est_max cat_price grademin grademax usonly novalue ringwin
gen targetobs = 0
save data/ps3_bidders, replace
restore


// Create data on target prices and append to knockout bid data
preserve
keep if nbidders == 2

drop if uniquelots == 0
sort auction
tab auction, gen(auc_)
replace bid = target_price
keep bid bidder auction auc_* est_min est_max cat_price grademin grademax usonly novalue ringwin
gen targetobs = 1

forval i=1/11{
	gen rm_`i' = 0
}

append using data/ps3_bidders

egen maxbid = max(bid) if targetobs != 1, by(auction)
gen ident = 0 if maxbid == bid
replace ident = 1 if maxbid != bid & targetobs != 1
replace ident = 2 if targetobs == 1
replace bidder = . if targetobs == 1
drop maxbid

save data/ps3_clean, replace
outsheet using data/ps3_clean.csv, comma nolabel replace

restore




*--------------------------------------------
* 		   			 Table 1
*-------------------------------------------- 

preserve

// collapse
duplicates drop house lot date, force
collapse (mean) target_mean = target_price (sd) target_sd = target_price ///
	(mean) bid_mean = mnbid (sd) bid_sd = mnbid  ///
	(mean) ringwin (sum) value ringvalue uniquelots, by(house)
	
// adjust percentages
gen percvalue = ringvalue/value*100
replace ringwin = ringwin*100
drop ringvalue value 

// order variables
local allvars1 target_mean target_sd bid_mean bid_sd ringwin percvalue uniquelots
order house `allvars1'

// Write to Latex
global Nval = _N
cap file close texfile
file open texfile using "${filepath}/tables/t1.tex", write replace
file write texfile "\begin{tabular}{lccccccc}" _n
file write texfile "\toprule" _n
file write texfile "& \multicolumn{2}{c}{Target Auction} & \multicolumn{2}{c}{Knockout Auction} & & & \\" _n 
file write texfile "\cmidrule(lr){2-3}\cmidrule(lr){4-5}" _n
file write texfile "Auction House & Mean & SD & Mean & SD & \% Lots won by Ring & \% Value Won & Total Lots \\ \midrule" _n
forval idx = 1/$Nval{
	local writer "`:di house[`idx']'"
	foreach var of varlist target_mean target_sd bid_mean bid_sd{
		local varval = trim("`: di %10.2f `var'[`idx']'")
		local writer "`writer' & `varval'"
		}
	foreach var of varlist ringwin percvalue uniquelots{
	local varval = trim("`: di %10.2g `var'[`idx']'")
		local writer "`writer' & `varval'"
	}
	file write texfile "`writer' \\" _n 
	}
file write texfile "\bottomrule" _n
file write texfile "\end{tabular}" _n
file close texfile
restore


*--------------------------------------------
* 		   			 Table 2
*-------------------------------------------- 

preserve

// collapse
duplicates drop house lot date, force
collapse (mean) target_mean = target_price (sd) target_sd = target_price ///
	(mean) bid_mean = medbid (sd) bid_sd = medbid  ///
	(mean) ringwin (sum) uniquelots, by(nbidders)
	
// order variables
replace ringwin = ringwin*100
local allvars2 target_mean target_sd bid_mean bid_sd ringwin uniquelots
order nbidders `allvars2'


// Write to Latex
global Nval = _N
cap file close texfile
file open texfile using "${filepath}/tables/t2.tex", write replace
file write texfile "\begin{tabular}{lcccccc}" _n
file write texfile "\toprule" _n
file write texfile "& \multicolumn{2}{c}{Target Auction} & \multicolumn{2}{c}{Knockout Auction} & & \\" _n 
file write texfile "\cmidrule(lr){2-3}\cmidrule(lr){4-5}" _n
file write texfile "\# Bidders & Mean & SD & Mean & SD & \% Lots won by Ring & Total Lots \\ \midrule" _n
forval idx = 1/$Nval{
	local writer "`:di nbidders[`idx']'"
	foreach var of varlist target_mean target_sd bid_mean bid_sd{
		local varval = trim("`: di %10.2f `var'[`idx']'")
		local writer "`writer' & `varval'"
		}
	foreach var of varlist ringwin uniquelots{
		local varval = trim("`: di %10.2g `var'[`idx']'")
		local writer "`writer' & `varval'"
		}
	file write texfile "`writer' \\" _n 
	}
file write texfile "\bottomrule" _n
file write texfile "\end{tabular}" _n
file close texfile
restore


*--------------------------------------------
* 		   			 Table 5
*-------------------------------------------- 

preserve


// collapse
*duplicates drop house lot date, force
collapse (mean) bidder_win (count) auction (mean) bidder_win2  ///
	(mean) sidepay payside (count) auction2, by(bidder)
	
// order variables
foreach var of varlist bidder_win bidder_win2 sidepay payside{
	replace `var' = `var'*100
}
local allvars5 bidder_win auction bidder_win2 sidepay payside auction2
order bidder `allvars5'

// Write to Latex
global Nval = _N
cap file close texfile
file open texfile using "${filepath}/tables/t5.tex", write replace
file write texfile "\begin{tabular}{b{1.5cm} b{1.5cm}b{2cm}b{1.5cm}b{2.3cm}b{2.3cm}b{2cm})}" _n
file write texfile "\toprule" _n
file write texfile "& \multicolumn{2}{c}{All Auctions $ (n \geq 1)$} & \multicolumn{4}{c}{Auctions w 2+ Ring Members $ (n \geq 2)$} \\" _n 
file write texfile "\cmidrule(lr){2-3}\cmidrule(lr){4-7}" _n
file write texfile "Ring Member & \% High KO Bid & \# Knockouts & \% High KO Bid & \% Receive Sidepayment & \% Pay Sidepayment & \# Knockouts \\ \midrule" _n
forval idx = 1/$Nval{
	local writer "`:di bidder[`idx']'"
	foreach var of varlist `allvars5'{
		local varval = trim("`: di %10.0f `var'[`idx']'")
		local writer "`writer' & `varval'"
		}
	file write texfile "`writer' \\" _n 
	}
file write texfile "\bottomrule" _n
file write texfile "\end{tabular}" _n
file close texfile
restore

*--------------------------------------------
* 		   			 Figure 1
*-------------------------------------------- 

preserve
gen netpay2 = netpay if target_price < 10000
label var bidder "Bidder"

graph bar (sum) netpayment netpay2, over(bidder) graphregion(color(white)) ///
	bar(1, bcolor(cranberry)) bar(2, bcolor(eltblue)) bargap(5) b1title("Bidder") ///
	ytitle("Net Sidepayment in Dollars") legend(size(small) label(1 "Net Sidepayment") ///
	label(2 "Net Sidepayment When Target Price < $10,000"))
graph export figures/fig1.png, replace
	
restore



