#####################################################
#         Metrics Problem Set 2 Question 2
#                     Tanya Rajan
# Description: This code runs lasso and ridge 
# regressions over 10,000 Monte Carlo draws 
# generated from the DGP specified in pset 1.
#####################################################

# setup
rm(list=ls())
library(data.table)
library(ggplot2)
library(knitr)
library(kableExtra)
library(MASS)
library(abind)
library(doParallel)
library(gsubfn)
library(pracma)
library(mvnfast)

# set working directory (change to run)
filepath <- "/Users/tanyarajan/Documents/git/y2/spring/metrics/pset2" 
setwd(filepath) 

# set seed
seeder = 123
set.seed(seeder, kind = "L'Ecuyer-CMRG")

# read in code file defining functions
source("code/pset2fns.R")

# setting cores for parallelization
cor<-floor(detectCores(all.tests=FALSE)*.5)
if (is.na(cor) | cor==0){cor <- 1}
registerDoParallel(cores=cor)





