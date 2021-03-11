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

// general cleaning
egen nbidders = count(bidder), by(lot house date)
egen totprofit = total(profit), by(lot house date)
egen highbid = max(bid), by(lot house date)
egen medbid = median(bid), by(lot house date)
gen ringwin = (target_price <= highbid)
egen uniquelots = tag(lot date house)
replace uniquelots = . if uniquelots == .
gen sidepay = (netpayment > 0)
gen payside = (netpayment < 0)
egen auction = group(lot date house)


// by house
egen h_target = mean(target_price) if uniquelots == 1, by(house)
egen h_target_sd = sd(target_price) if uniquelots == 1, by(house)
egen h_knock = mean(highbid) if uniquelots == 1, by(house)
egen h_knock_sd = sd(highbid) if uniquelots == 1, by(house)
egen h_percwin = mean(ringwin) if uniquelots == 1, by(house)
replace h_percwin = h_percwin*100
gen h_value = (est_min + est_max)/2 
egen h_totvalue = total(h_value) if uniquelots == 1, by(house) missing
egen h_ringvalue = total(h_value) if uniquelots == 1 &  ringwin == 1, by(house) missing
gen h_percvalue = h_ringvalue/h_totval * 100
egen h_totlots = total(uniquelots), by(house) missing


// by number of bidders
*egen uniquelots = tag(lot date house)
egen b_target = mean(target_price) if uniquelots == 1, by(nbidders)
egen b_target_sd = sd(target_price) if uniquelots == 1, by(nbidders)
egen b_knock = mean(medbid) if uniquelots == 1, by(nbidders)
egen b_knock_sd = sd(highbid) if uniquelots == 1, by(nbidders)
egen b_percwin = mean(ringwin) if uniquelots == 1, by(nbidders)
replace b_percwin = b_percwin*100
egen b_totlots = total(uniquelots), by(nbidders) missing


// by ring member
gen m_wins = (rank == 1)
egen m_percwins = mean(m_wins), by(bidder)
replace m_percwins = m_percwins*100 
egen m_totlots = total(uniquelots), by(bidder) missing


// by ring member if auction has more than 2 bidders
gen m2_wins = (rank == 1 & nbidders >= 2)
egen m2_percwins = mean(m2_wins), by(bidder)
replace m2_percwins = m2_percwins*100 
egen m2_totlots = total(uniquelots) if nbidders >=2, by(bidder) missing
egen m2_side = mean(sidepay) if nbidders >=2, by(bidder)
replace m2_side = m2_side*100
egen m2_pay = mean(payside) if nbidders >=2, by(bidder)
replace m2_pay = m2_pay*100


foreach var of varlist h_* b_* m_* m2_*{
	replace `var' = . if uniquelots == .
}

// all manipulations
save data/allvars, replace


*** Keep only data needed for analysis ***

// Data on knockout bids
preserve
drop b_* h_* m2_* 
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
drop b_* h_* m2_* 
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

// Keep only one obs per house
*drop if h_ringvalue == .
duplicates drop house, force
sort house
keep house h_*
drop h_ringvalue h_value h_totvalue


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
	foreach var of varlist h_*{
		local varval = trim("`: di %10.3g `var'[`idx']'")
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

// Keep only one obs per house
drop if b_percwin == .
duplicates drop nbidders, force
sort nbidders
keep nbidders b_*

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
	foreach var of varlist b_*{
		local varval = trim("`: di %10.3g `var'[`idx']'")
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

// Keep only one obs per house
drop if uniquelots == .
duplicates drop bidder, force
sort bidder
keep bidder m_* m2_*
drop m_wins m2_wins
order bidder m_percwins m_totlots m2_percwins m2_side m2_pay m2_totlots

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
	foreach var of varlist m*{
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

graph bar netpayment netpay2, over(bidder) graphregion(color(white)) ///
	bar(1, bcolor(cranberry)) bar(2, bcolor(eltblue)) bargap(5) b1title("Bidder") ///
	ytitle("Net Sidepayment in Dollars") legend(size(small) label(1 "Net Sidepayment") ///
	label(2 "Net Sidepayment When Target Price < $10,000"))
graph export figures/fig1.png, replace
	
restore



