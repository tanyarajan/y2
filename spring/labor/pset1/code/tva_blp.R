#####################################################
#               Labor Problem Set 1
#                     Tanya Rajan
# Description: This code runs the BLP part of the 
# TVA estimation
#####################################################

# setup
rm(list=ls())
library(data.table)
library(ggplot2)
library(knitr)
library(kableExtra)

# set working directory (change to run)
filepath <- "/Users/tanyarajan/Documents/git/Education/pset1" 
setwd(filepath) 

# read in data
dfbig <- as.data.table(read.csv(file = "data/va_clean.csv"))

# function to evaluate arguments in text
callit <- function(...){
  arg <- parse(text=paste0(..., collapse=""))
  eval(parse(text=arg))
}

#####################################################
#         Calculating Variance Components           #
#####################################################
# degrees of freedom correction from appendix
N <- uniqueN(dfbig$id_student)
C <- uniqueN(dfbig$class)
K <- 60 # from Stata: number of covariates used
dfcorr <- (N - 1)/(N - K - C + 1)

# calculating total variance of A
s_A <- var(dfbig$A, na.rm=T)*dfcorr

# calculating variance of deviations from class mean
s_ep <- var(dfbig$A - dfbig$class_A, na.rm=T)*dfcorr

# calculating (s_A0 + s_theta)
var_class <- s_A - s_ep


#####################################################
#               Creating BLP matrices               #
#####################################################
# read in the collapsed data
df <- as.data.table(read.csv(file = "data/va_collapseclean.csv"))

# create selection matrix based on how many available lags
selcols = c("A", "A_l1", "A_l2", "A_l3", "A_l4", "A_l5", "A_l6")
selection_idx = !is.na(df[,..selcols])

# all possible values of gamma (for all lags)
G <- list()
for (i in 1:length(selcols)){
  covs[i] <- var(df$A, callit("df$A_l", i), na.rm=T)
}
G <- unlist(covs)

# creating components of Sigma common across all teachers
Sig_init = var(df[,..selcols], na.rm = T)
Sig_init2 = Sig_init - diag(diag(Sig_init)) + diag(rep(var_class, length(G)))

# function to calculate teacher value add
tva <- function(idx){
  # relevant selection vectors
  selmat <- diag(selection_idx[idx,])
  Avec <-t(df[idx,..selcols])
  Avec[is.na(Avec)] <- 0
  
  # relevant number of students per teacher-year
  nct <- df[idx, ns]
  weights <- diag(rep(s_ep/nct, length(selection_idx[idx,])))
  
  # constructing sigma/gamma
  gamma <- selmat%*%G
  Sigma <- selmat%*%Sig_init2%*%selmat + weights
  
  # coeff
  psi <- solve(Sigma)%*%gamma
  psi[is.na(psi)] <- 0
  
  return(t(psi)%*%Avec)
}

# calculating tva ests for all teacher-years
tva_ests <- sapply(rep(1:dim(df)[1]), tva)
df$tva_ests <- tva_ests

# normalizing and plotting
tva_ests <- (tva_ests - mean(tva_ests))/sd(tva_ests)
plot(density(tva_ests, bw=2))