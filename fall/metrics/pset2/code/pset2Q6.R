#####################################################
#               Metrics Problem Set 2
#                      Question 6
#                     Tanya Rajan
# Description: This code runs Monte Carlo simulations
# for various estimators given the DGP stated in the
# question.
#####################################################

# setup
rm(list=ls())
require(data.table)
require(ggplot2)
require(pracma)
require(doParallel)
require(MASS)
require(SteinIV)
require(knitr)
require(tidyr)
require(kableExtra)

# set working directory (change to run)
filepath <- "/Users/tanyarajan/Documents/git/psets/metrics_fall2020/pset2" 
setwd(filepath) 

# call in code file with function definitions
source("code/pset2fns.R")


############ DGP and Simulation Parameters ##############
set.seed(1234)


simulate <- function(N,D,mumat, varmat){
  out <- list()
  
  # drawing U,V jointly
  mcmat <- lapply(1:D, function(x) mvrnorm(N, mumat, varmat))
  
  # drawing Z
  zmat <- lapply(1:D, function(x) matrix(rnorm(N*20, 0, 1),ncol = 20) )
  out$z <- zmat
  
  # calculating rest of data
  zprep <- lapply(zmat, function(x) x[,1]*0.3) 
  vprep <- lapply(mcmat, function(x) x[,2])
  xmat <- Map("+", zprep, vprep)
  out$x <- xmat
  
  uprep <- lapply(mcmat, function(x) x[,1])
  ymat <- Map("+", xmat, uprep)
  out$y <- ymat
  
  # return data
  return(out)
  
}


# mean and variance parameters
mu <- matrix(c(0,0), ncol=1)
var <- matrix(c(.25, .20, .20, .25), ncol=2)


########### Running the simulations ##############
draws <- 100 # set low to run code faster
outarray <- array(NA,dim=c(draws,5,4))
outcovs <- array(NA, dim=c(draws,5,4))
c = 1
for (nn in c(100, 200, 400, 800)){
  
  # simulate data
  sim <- simulate(nn,draws,mu,var)
  ymat <- sim$y
  xmat <- sim$x
  zmat <- sim$z
  
  for (t in c("ols", "tsls", "jack", "tslsmult", "jackmult") ){
    callit("est_",t, "<<-list()")
    callit("cov_",t, "<<-list()")
  }
  
  for (d in 1:draws){
    out_ols <- reg(ymat[[d]], xmat[[d]])
    out_tsls <- iv(ymat[[d]], xmat[[d]], zmat[[d]][,1])
    out_jack <- iv(ymat[[d]], xmat[[d]], zmat[[d]][,1], est="jack")
    out_tslsmult <- iv(ymat[[d]], xmat[[d]], zmat[[d]][,1:20])
    out_jackmult <- iv(ymat[[d]], xmat[[d]], zmat[[d]][,1:20], est="jack")
    
    for (t in c("ols", "tsls", "jack", "tslsmult", "jackmult") ){
      callit("est_",t,"[d] <<- out_",t,"$b[2]")
      callit("cov_",t,"[d] <<- out_",t,"$cov")
    }
  
  }
  
  # store results
  outarray[,,c] <- cbind(unlist(est_ols), unlist(est_tsls), unlist(est_jack),
                        unlist(est_tslsmult), unlist(est_jackmult))
  outcovs[,,c] <- cbind(unlist(cov_ols), unlist(cov_tsls), unlist(cov_jack),
                         unlist(cov_tslsmult), unlist(cov_jackmult))
  c = c+1
  
}


########### Outputting to Latex Table ##############
outmeds <- apply(outarray, c(2,3), function(x) quantile(x,.5))
outsd <- apply(outarray, c(2,3), function(x) sd(x))
outbias <- apply(outarray, c(2,3), function(x) mean(x-1))
outcoverage <- apply(outcovs, c(2,3), function(x) mean(x))
outtable <- c(NA, NA, NA, NA)
for (m in 1:5){
  outtable <- rbind(outtable, round(outmeds[m,],4))
  outtable <- rbind(outtable, round(outbias[m,],4))
  outtable <- rbind(outtable, round(outsd[m,],4))
  outtable <- rbind(outtable, round(outcoverage[m,],4))
}
outtable<-outtable[-1,]
names <- c("Median", "Bias","SD", "Coverage", 
           "Median",  "Bias","SD", "Coverage",
           "Median", "Bias","SD", "Coverage",
           "Median", "Bias","SD", "Coverage", 
           "Median", "Bias","SD", "Coverage")
outtable <- cbind(names, outtable)
colnames(outtable)<-c(" ", "N=100","N=200","N=400","N=800")

kable(outtable, format="latex", booktabs=TRUE) %>%
  add_header_above(c(" "=2, "Sample Sizes" = 4)) %>%
  pack_rows("OLS, No Instruments",1,4) %>%
  pack_rows("TSLS, 1 Instrument",5,8) %>%
  pack_rows("Jackknife TSLS, 1 Instrument",9,12) %>%
  pack_rows("TSLS, Many Instruments",13,16) %>%
  pack_rows("Jackknife TSLS, Many Instruments",17,20) %>%
  write(file="tables/tr_MCmeans.tex")




