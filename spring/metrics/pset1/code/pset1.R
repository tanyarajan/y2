#####################################################
#               Metrics Problem Set 1
#                     Tanya Rajan
# Description: This code runs OLS with varying
# numbers of covariates across Monte Carlo samples
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
library(parallel)
library(future.apply)
library(pracma)

# set working directory (change to run)
filepath <- "/Users/tanyarajan/Documents/git/y2/spring/metrics/pset1" 
setwd(filepath) 

# set seed
seeder = 123
set.seed(seeder)

# read in code file defining functions
source("code/pset1fns.R")

# setting cores for parallelization
cor<-floor(detectCores(all.tests=FALSE)*.7)
if (is.na(cor) | cor==0){cor <- 1}
registerDoParallel(cores=cor)
plan(multisession, workers=cor)

#####################################################
#                   Data Generation                 #
#####################################################

drawdata <- function(N, P, M, vcv="ident", rho=NULL){
  # mean 
  mu = rep(0, P+1)
  
  # different variance matrix (either identity or autocovarying)
  if (vcv == "ident"){sigma=eye(P+1)}
  else{
    sigma1 = eye(P+1)*(rho - 1)
    sigma2 = matrix(rho, nrow=P+1, ncol=P+1)
    sigma = sigma2 - sigma1
      }
  
  # drawing x, epsilon, and y
  xep = abind(future_lapply(1:M, function(x){mvrnorm(N, mu, sigma)}, future.seed=seeder), along=3)
  y = xep[,1,] - xep[,2,] + xep[,P+1,]
  x = xep[,-(P+1),]
  return(list(x=x, y=y))
}


#####################################################
#     Functions to Compute Various Objects          #
#####################################################

# Calculating bhat_1 
preg_yx <- function(data,p,M){
  looper <- foreach (j=1:M) %dopar% {
    reg(data$y[,j], data$x[,1:p,j], const="none")$b[1]
  }
  return(unlist(looper))
}

# Calculating xhat_1 
preg_xx <- function(data,p,M){
  looper <- foreach (j=1:M) %dopar% {
    mean((data$x[,2:p,j]%*%reg(data$x[,1,j], data$x[,2:p,j], const="none")$b)^2)
  }
  return(unlist(looper))
}

# Calculating eigenvalues of X'X
eigenx <- function(data, p, M, ret="min"){
  looper <- foreach (j = 1:M) %dopar% {
    return(sort(eigen(t(data$x[,2:p,j])%*%data$x[,2:p,j], only.values=T)$values))
  }
  loopres <- abind(looper, along=2)
  if (ret=="min"){return(mean(loopres[1,]))}
  else{return(rowMeans(loopres))}
}

#####################################################
#                Questions b - d                    #
#####################################################
setM <- 10000
setP <- 90
setN <- 100
plist <- c(1,5,10,50,85,90)

# Drawing the sample for part b-e
sample1 <- drawdata(setN, setP, setM)

# Calculating means and variances for part c
bhats1 <- lapply(plist, function(p){preg_yx(sample1, p, setM)})
bmean1 <- sapply(bhats1, mean)
bse1 <- sapply(bhats1, var)

# Calculating means and variances for part d
plist <- c(5,10,50,85,90)
xhats1 <- lapply(plist, function(p){preg_xx(sample1, p, setM)})
xmeans1 <- sapply(xhats1, mean)

# Calculating eigenvalues for part d
eigens1 <- lapply(plist, function(p){eigenx(sample1, p, setM)})
eigmean1 <- sapply(eigens1, mean)

#####################################################
#                  Questions e - f                  #
#####################################################

Nlist <- c(100, 200, 500, 1000)
setM2 <- 1000

# function to do all parts of questions e and f
output_e <- function(M, Nvals, P, rho, pdep=F, fn=NULL){
  if (pdep==T){Pvals = sapply(Nvals, fn)}
  else{Pvals = c(P,P,P,P)}
  samples_e <- lapply(1:4, function(n){drawdata(Nvals[n], Pvals[n], M, vcv="auto", rho=rho)})
  
  # average and variance of bhat_1
  bhat <- lapply(1:4, function(n){preg_yx(samples_e[[n]], Pvals[n], M)})
  bmean <- sapply(bhat, mean)
  bse <- sapply(bhat, var)
  
  # lowest eigenvalues
  eigen <- lapply(1:4, function(n){eigenx(samples_e[[n]], Pvals[n],M)})
  eigmean <- sapply(eigen, mean)
  
  # plot averages of the ordered eigenvalues
  eigenplot <- eigenx(samples_e[[4]], Pvals[n], M, ret="all")
  
  return(list(bmn=bmean, bse=bse, eigmean=eigmean, eigplot=eigenplot))
  
}


out_rho0 <- output_e(setM, Nlist, setP, 0)
out_rho5 <- output_e(setM, Nlist, setP, .5)
out_rho9 <- output_e(setM, Nlist, setP, .9)


#####################################################
#         Outputting Tables and Graphs              #
#####################################################
# Table 1
table1 <- rbind(bmean1, bse1, c(NA, xmeans1), c(NA, eigmean1))
rownames(table1) <- c("mean($\\hat{\\beta}$)", "var($\\hat{\\beta}$)", 
                      "mean($\\frac{1}{N}\\sum \\hat{X}_{i1}(\\tilde{p})^2$)",
                      "mean($\\underbar{s}_k$)")
colnames(table1) <- c("$\\tilde{p}=1$", "$\\tilde{p}=5$", "$\\tilde{p}=10$",
                      "$\\tilde{p}=50$", "$\\tilde{p}=85$", "$\\tilde{p}=90$")
table1 <- round(table1, 4)

# Write to Latex
kable(table1, format="latex", booktabs=TRUE, escape=FALSE) %>%
  write(file="tables/pset1_table1.tex")

# Table 2
table2 <- rbind(out_rho0$bmn, out_rho0$bse, out_rho0$eigmean,
                out_rho5$bmn, out_rho5$bse, out_rho5$eigmean,
                out_rho9$bmn, out_rho9$bse, out_rho9$eigmean)
colnames(table2) <- c("N=100", "N=200", "N=500", "N=1000")
rownames(table2) <- rep(c("mean($\\hat{\\beta}$)", "var($\\hat{\\beta}$)", 
                          "mean($\\underbar{s}_k$)"), 3)
table2 <- round(table2, 4)

# Write to Latex
kable(table2, format="latex", booktabs=T, escape=FALSE) %>%
  add_header_above(c(" "=1, "Sample Sizes" = 4)) %>%
  pack_rows("$ \\rho = 0 $  ", 1, 3, escape=F) %>%
  pack_rows("$ \\rho = 0.5 $", 4, 6, escape=F) %>%
  pack_rows("$ \\rho = 0.9 $", 7, 9, escape=F) %>% 
  write(file="tables/pset1_table2.tex")


# Table 3
table3 <- rbind(outg_rho0$bmn, outg_rho0$bse, outg_rho0$eigmean,
                outg_rho5$bmn, outg_rho5$bse, outg_rho5$eigmean,
                outg_rho9$bmn, outg_rho9$bse, outg_rho9$eigmean)
colnames(table3) <- c("N=100", "N=200", "N=500", "N=1000")
rownames(table3) <- rep(c("mean($\\hat{\\beta}$)", "var($\\hat{\\beta}$)", 
                          "mean($\\underbar{s}_k$)"), 3)
table3 <- round(table3, 4)

# Write to Latex
kable(table3, format="latex", booktabs=T, escape=FALSE) %>%
  add_header_above(c(" "=1, "Sample Sizes" = 4)) %>%
  pack_rows("$ \\rho = 0 $  ", 1, 3, escape=F) %>%
  pack_rows("$ \\rho = 0.5 $", 4, 6, escape=F) %>%
  pack_rows("$ \\rho = 0.9 $", 7, 9, escape=F) %>% 
  write(file="tables/pset1_table3.tex")



# Figure 1
png(file="figures/pset1_fig1.png", width=600, height=600)
par(mfrow=c(2,2))
plot(seq(1:length(out_rho0$eigplot)), rev(out_rho0$eigplot), type="l", 
     main=expression(paste(rho, " = 0")), xlab="Eigenvalue Index", ylab="Eigenvalue")
plot(seq(1:length(out_rho5$eigplot)), rev(out_rho5$eigplot), type="l", 
           main=expression(paste(rho, " = 0.5")), xlab="Eigenvalue Index", ylab="Eigenvalue")
plot(seq(1:length(out_rho9$eigplot)), rev(out_rho9$eigplot), type="l", 
           main=expression(paste(rho, " = 0.9")), xlab="Eigenvalue Index", ylab="Eigenvalue")
dev.off()


# Figure 1b zoomed in
eigendata <- as.data.table(cbind(seq(1:length(out_rho0$eigplot))[-1], 
                                 rev(out_rho0$eigplot)[-1], 
                                 rev(out_rho5$eigplot)[-1], 
                                 rev(out_rho9$eigplot)[-1]))
names(eigendata) <- c("idx", "rho = 0", "rho = 0.5", "rho = 0.9")
eigenlong <- melt(eigendata, id=c("idx"))
g2<-ggplot(data=eigenlong, aes(x=idx, y=value, colour=variable)) + geom_line() + theme_bw()
ggsave(g2, width = 7, height = 5, file = "figures/pset1_fig1b.png")




