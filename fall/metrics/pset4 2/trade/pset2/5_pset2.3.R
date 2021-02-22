#####################################################
#     International Macro and Trade Problem Set 2
#                     Tanya Rajan
# Description: This file uses data associated with Head
# and Mayer (2014) to run the fixed-effects 
# log-linear version of the gravity regression on all
# years from 1948-2006.
#####################################################

# setup
rm(list=ls())
require(data.table)
require(ggplot2)
require(foreign)
require(fixest)
require(haven)

# set working directory
filepath <- "/Users/tanyarajan/Documents/git/psets/trade/pset2" # change to run
tables <- paste0(filepath,"/tables")
setwd(filepath) 

#################### Estimation ###################
# load data
data <- as.data.table(read_dta("data/clean.dta"))
subset <- data[flow!=0,]
setnames(subset, "lndist", "LogDistance")
setnames(subset, "lnflow", "LogFlow")
setnames(subset, "comlang_off", "CommonLanguage")
setnames(subset, "contig", "Contiguous")

# using the fixest package
time<-system.time(
output <- feols(data=subset, LogFlow ~ LogDistance +
                  Contiguous + CommonLanguage | expyear +impyear )
)

# write to latex using heteroskedasticity-robust errors
etable(output, se="white", drop.section=list("fixef", "stats"))




