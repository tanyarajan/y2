#####################################################
#               Metrics Problem Set 4
#                      Question 2
#                     Tanya Rajan
# Description: This file runs code for Question 2
# on the problem set.
#####################################################

# setup
rm(list=ls())
require(data.table)
require(ggplot2)
require(pracma)
require(nloptr)
require(optimx)
require(knitr)
require(kableExtra)

# set working directory (change to run)
filepath <- "/Users/tanyarajan/Documents/git/psets/metrics_fall2020/pset4" 
setwd(filepath) 

# read in code file defining functions
source("code/pset4fns.R")

# null
nullhyp <- sin(1)


#####################################################
#           Data Preparation Functions              #
#####################################################
# Drawing errors for all Monte carlo iterations
drawDGP <- function(n, r, m){
  set.seed(1234)
  err.E  <- lapply(1:m, function(x){sample.int(4, size=n, prob = rep(1/4, 4), replace=T) + 1})
  err.V  <- lapply(1:m, function(x){matrix(rnorm(n*5), nrow=n, ncol=5)})
  err.ep <- lapply(1:m, function(x){matrix(rnorm(n*5), nrow=n, ncol=5)})
  err.U <- list()
  for (k in 1:m){
    U <- matrix(err.ep[[k]][,1], ncol=1)
    for (i in 2:5){
      idx = i - 1
      U <- cbind(U, r*U[,idx] + err.ep[[k]][,i])
    }
    err.U[[k]] <- U
  }
  return(list(E = err.E, V=err.V, U = err.U))
}


# Function to generate "observed" data given parameters
gendata <- function(n, th, E, V, U){
    
    # Calculating Y(0), Y(1)
    Y0 <- -0.2 + 0.5*E + U
    Y1 <- sapply(1:5, function(i){-0.2 + 0.5*E + sin(i - th*E) + U[,i] + V[,i]  })
    colnames(Y0) <- c(1,2,3,4,5)
    colnames(Y1) <- c(1,2,3,4,5)
      
    # Creating long dataset
    Y0L <- melt(as.data.table(cbind(1:n, Y0)), id.vars="V1")
    Y1L <- melt(as.data.table(cbind(1:n, Y1)), id.vars="V1")
    gendata <- as.data.table(cbind(Y0L, Y1L[,c(-1,-2)], rep(E,5)))
    names(gendata) <- c("person", "time", "Y0", "Y1", "E")

    # Outcomes
    gendata[,D:= as.numeric(time) >= E]
    gendata[,Y:=Y0*(1-D) + Y1*D]
    
    # Cohort FE
    gendata[,`:=`(E2 = E==2, E3 = E == 3, E4 = E == 4, E5 = E == 5)]
    
    # Time FE
    gendata[,`:=`(T1 = time==1,T2 = time==2,T3 = time==3,T4 = time==4, T5 = time==5)]
    
    # Relative Time FE
    rellist <- c(-3, -2, 0, 1, 2, 3)
    gendata[,reltime:=as.numeric(time)-as.numeric(E)]
    idx = 1
    for (r in rellist){
      gendata[,paste0("rel",idx):= reltime == r]
      idx = idx+1
    }
    return(gendata)
    
}

# Function recovering regression parameters / ATE
regdata <- function(gendata, obj="reg", setype="", reps=NULL){
    
    # Running the regression on relative time indicators
    selx <- c("E3", "E4", "E5", "T2", "T3", "T4", "T5", 
              "rel1", "rel2", "rel3", "rel4", "rel5", "rel6")
    yreg <- as.matrix(gendata$Y, ncol=1)
    xreg <- as.matrix(subset(gendata, select=selx))
    if (setype == "cluster"){cluster <- as.matrix(gendata$person, ncol=1)}
    if (setype != "cluster"){cluster <- NULL}
    regout <- reg(yreg,xreg, se=setype, cluster)
    regcoeff <- regout$b
    regse <- regout$se
    rt.coeffs <- rev(rev(regcoeff)[1:6])
    rt.se <- rev(rev(regse)[1:6])
    
    # regression estimates
    if (obj == "reg"){
      return(rt.coeffs)
    }
    
    # ATE estimates and actual ATE
    if (obj == "ate"){
      gendata[,diff:=Y1 - Y0]
      ate <- gendata[E==2 & time==3, mean(diff)]
      term1 <- gendata[E == 2 & time == 3, mean(Y)]
      baseline <- gendata[E == 2 & time == 1, mean(Y)]
      ct <- gendata[E==4 & time==3, mean(Y)] - gendata[E==4 & time==1, mean(Y)]
      ateest <- term1 - baseline - ct
      return(c(ate, ateest))
    }
    
    # Hypothesis test on first coefficient
    if (obj == "test"){
      stat <- (rt.coeffs[4] - nullhyp)/(rt.se[4])
      return((abs(stat) > qnorm(.975)))
    }
}

# Plotting function
mc.plotter <- function(small, large){
  for (i in c("small", "large")){
    assign(paste0("c.", i),  matrix(unlist(get(i)), ncol = 6, byrow=T))
    assign(paste0("mn.",i), colMeans(get(paste0("c.",i))))
    assign(paste0("up.",i), apply(get(paste0("c.",i)), 2, function(x){quantile(x, c(.975))}))
    assign(paste0("low.",i), apply(get(paste0("c.",i)), 2, function(x){quantile(x, c(.025))}))
  }
  reltime <- c(-3, -2, 0, 1, 2, 3)
  plotter <- as.data.table(cbind(reltime, mn.small, mn.large, up.small, up.large, low.small, low.large))
  g <- ggplot(data = plotter, aes(x=reltime)) + 
    geom_ribbon(aes(ymin=low.small, ymax=up.small), fill = "coral2", alpha=.4) + 
    geom_ribbon(aes(ymin=low.large, ymax=up.large), fill = "cornflowerblue", alpha=.4) + 
    geom_line(aes(y=mn.small, color="N=1000"))  + 
    geom_line(aes(y=mn.large, color="N=10000")) +
    scale_color_manual(values = c("N=1000"="coral2", "N=10000"="cornflowerblue")) + 
    theme_bw() + ylab("Coefficient") + xlab("Relative Time Indicator") + labs(color = 'Number of Obs.')
  return(g)
}

#####################################################
#              Running Regressions                  #
#####################################################

# Monte Carlo parameters
M <- 500

# Drawing all variations
draw1<-drawDGP(n=1000, r=0.5, m=M)
draw2<-drawDGP(n=10000, r=0.5, m=M)

# generate data
data11 <- lapply(1:M, function(k){gendata(n=1000, th=-2, draw1$E[[k]], draw1$V[[k]],draw1$U[[k]])} )
data12 <- lapply(1:M, function(k){gendata(n=10000, th=-2, draw2$E[[k]], draw2$V[[k]],draw2$U[[k]])})
data21 <- lapply(1:M, function(k){gendata(n=1000, th=0, draw1$E[[k]], draw1$V[[k]],draw1$U[[k]])})
data22 <- lapply(1:M, function(k){gendata(n=10000, th=0, draw2$E[[k]], draw2$V[[k]],draw2$U[[k]])})
data31 <- lapply(1:M, function(k){gendata(n=1000, th=1, draw1$E[[k]], draw1$V[[k]],draw1$U[[k]])})
data32 <- lapply(1:M, function(k){gendata(n=10000, th=1, draw2$E[[k]], draw2$V[[k]],draw2$U[[k]])})

# Getting regression coefficients
reg11 <- lapply(1:M, function(k){regdata(data11[[k]])})
reg12 <- lapply(1:M, function(k){regdata(data12[[k]])})
reg21 <- lapply(1:M, function(k){regdata(data21[[k]])})
reg22 <- lapply(1:M, function(k){regdata(data22[[k]])})
reg31 <- lapply(1:M, function(k){regdata(data31[[k]])})
reg32 <- lapply(1:M, function(k){regdata(data32[[k]])})

# Plotting
theta_2 <- mc.plotter(reg11, reg12)
ggsave(theta_2, width = 7, height = 5, file = "figures/tr_theta_2.png")
theta0 <- mc.plotter(reg21, reg22)
ggsave(theta0, width = 7, height = 5, file = "figures/tr_theta0.png")
theta1 <- mc.plotter(reg31, reg32)
ggsave(theta1, width = 7, height = 5, file = "figures/tr_theta1.png")

# ATEs
ate11 <- lapply(1:M, function(k){regdata(data11[[k]], obj="ate")})
ate12 <- lapply(1:M, function(k){regdata(data12[[k]], obj="ate")})
ate21 <- lapply(1:M, function(k){regdata(data21[[k]], obj="ate")})
ate22 <- lapply(1:M, function(k){regdata(data22[[k]], obj="ate")})
ate31 <- lapply(1:M, function(k){regdata(data31[[k]], obj="ate")})
ate32 <- lapply(1:M, function(k){regdata(data32[[k]], obj="ate")})

# Writing to Latex
atetable <- matrix(NA, nrow=3, ncol=2)
atetable[1,]<-colMeans(matrix(unlist(ate12), ncol = 2, byrow=T))
atetable[2,]<-colMeans(matrix(unlist(ate22), ncol = 2, byrow=T))
atetable[3,]<-colMeans(matrix(unlist(ate32), ncol = 2, byrow=T))
colnames(atetable) <- c("Actual", "Estimated")
rownames(atetable) <- c("$\\theta = -2$", "$\\theta = 0$", "$\\theta = 1$")

kable(atetable, digits=4, format="latex", booktabs=TRUE, escape=F) %>%
  write(file="tables/tr_ates.tex")



#####################################################
#                Part E: Testing                    #
#####################################################
# DGP parameters
N <- rep(c(40, 50, 200), each = 3)
rho <- rep(c(0, .5, 1), times = 3)
theta <- 1
M <- 500
B <- 100 # bootstrap reps (reduced here to run code faster)

# Empirical rejection rates
tester <- function(n, r, m, th, setype="", reps=NULL){
  draw <-  drawDGP(n=n, r=r, m=M)
  data <-  lapply(1:m, function(k){gendata(n=n, th=th, draw$E[[k]], draw$V[[k]],draw$U[[k]])} )
  if (setype != "wild"){
    test <-  lapply(1:m, function(i){regdata(data[[i]], obj="test", setype=setype)})
  }
  if (setype == "wild"){
    test <- lapply(1:m, function(i){wildboot(data[[i]], reps, true=nullhyp)})
  }
  rej  <-  mean(as.numeric(unlist(test)))
  return(rej)
}


# Calculating rejection rates by SE type
rej_se1 <- mapply(tester, N, rho, m=M, th=theta, setype="")
rej_se2 <- mapply(tester, N, rho, m=M, th=theta, setype="HC1")
rej_se3 <- mapply(tester, N, rho, m=M, th=theta, setype="cluster")
rej_se4 <- mapply(tester, N, rho, m=100, th=theta, setype="wild", reps=50) 

# Writing to Latex
rejtable <- rbind(matrix(rej_se1, nrow=3), matrix(rej_se2, nrow=3), 
      matrix(rej_se3, nrow=3), matrix(rej_se4, nrow=3))
colnames(rejtable) <- c("N=40", "N=50", "N=200")
rownames(rejtable) <- rep(c("$\\rho = 0$", "$\\rho = 0.5$", "$\\rho = 1$"), 4)

kable(rejtable, format="latex", booktabs=T, escape=F) %>%
  add_header_above(c(" "=1, "Sample Sizes" = 3)) %>%
  pack_rows("Homoskedastic SEs", 1, 3) %>%
  pack_rows("Robust SEs", 4, 6) %>%
  pack_rows("Clustered Robust SEs", 7, 9) %>%
  pack_rows("Clustered Wild Bootstrap", 10, 12)%>%
  write(file="tables/tr_rejrates.tex")







