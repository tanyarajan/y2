#####################################################
#               Metrics Problem Set 1
#                      Question 4
#                     Tanya Rajan
# Description: This code plots the results of a Monte
# Carlo simulation exercise for various nonparametric
# estimators.
#####################################################

# setup
rm(list=ls())
require(data.table)
require(ggplot2)
require(reshape2)
require(pracma)
require(sfsmisc)
require(OneR)
require(Hmisc)
require(doParallel)
require(gridExtra)

# set working directory (change to run)
filepath <- "/Users/tanyarajan/Documents/git/psets/metrics_fall2020/pset1" 
setwd(filepath) 

############## Simulation Parameters ##############

# set model parameters
meanU <- 0
sdU <- 0.3
upperX <- 2
lowerX <- -2

# set Monte Carlo parameters
set.seed(1234)
N <- 500 # number of observations per draw
D <- 100 # number of draws


# drawing observations
U <- matrix(rnorm(N*D, meanU, sdU), N, D)
X <- matrix(runif(N, lowerX, upperX), N, D) # X remains fixed across draws
Y <- apply(X, c(1,2), function(x) sin(2*x)) + 
  apply(X, c(1,2), function(x) 2*exp(-16*(x^2))) + U

# parallelization requirements (setting # of computer cores to use)
cor<-floor(detectCores(all.tests=FALSE)*.7)
if (is.na(cor) | cor==0){cor <- 1}
registerDoParallel(cores=cor)

# fonts
myFont <- "Palatino"


############# Functions for Each Method ############

# setting up kernels
kernel <- function(kern, arg){
  if (kern=="unif"){
    return((1/2)*(abs(arg) <= 1))
  }
  else if (kern == "gauss"){
    return(1/(sqrt(2*pi))*exp(-(arg^2)/2))
  }
  else if (kern == "epan"){
    return((3/4)*(1 - arg^2)*(abs(arg) <= 1))
  }
}

### local constant (kernel) regression ###
lcr <- function(h, y, x, kern){
  outres <- rep(NA, length(x))
  lcr_dta <- as.data.table(cbind(y,x))
  
  # for x within bandwidth, take average
  for (i in 1:length(x) ){
    val = x[i]
    lcr_dta[, weight:=kernel(kern, ((x - val)/h))]
    denom <- sum(lcr_dta$weight)
    num <- sum(lcr_dta$weight*lcr_dta$y)
    outres[i] <- num/denom
  }
  return(outres) 
}

### local linear regression ##
llr <- function(h, y, x, kern){
  outres <- rep(NA, length(x))
  llr_dta <- as.data.table(cbind(1,y,x))
  const <- lcr(h,y,x,kern)
  xw <- lcr(h,x,x,kern) # LCR but with x everywhere instead of y
  
  # for x within bandwidth, apply LLR estimator
  for (i in 1:length(x)){
    val = x[i]
    iter <- cbind(llr_dta$V1, llr_dta$x)
    llr_dta[, delta := x[i] - xw]
    llr_dta[, weight:=kernel(kern, ((x-val)/h) ) ]
    denom <- sum(llr_dta$weight*(llr_dta$delta)^2)
    num <- sum(llr_dta$weight*llr_dta$delta*llr_dta$y)
    outres[i] <- (num/denom)
  }
  outres2 <- const + (x - xw)*outres
  return(outres2)
}


### sieve approximation with standard polynomial basis ###
sieve <- function(K,y,x){
  sieve_dta <- as.data.table(cbind(y,1))
  
  # set up sieve polynomial and regress
  for (i in 1:K){
    sieve_dta[,paste0("x",i):=x^i]
  }
  Y<- as.matrix(y,ncol=1)
  X<- as.matrix(sieve_dta[,y:=NULL])
  ols<- solve(t(X)%*%X, tol = 1e-20)%*%(t(X)%*%Y)
  out <- X%*%ols
  return(out)
}

### nearest neighbors estimator ###
nn <- function(K,y,x,match){
  outres <- rep(NA, length(x))
  nn_dta <- as.data.table(cbind(y,x))
  if (match=="maha"){
    w<-var(x) # calculate variance for mahalanobis
  }
   else{
     w<-1
   }
  
  # find nearest neighbors
  for (i in 1:length(x)){
    nn_dta[, dist:=((x - x[i])^2)/w]
    mn<-head(nn_dta[order(dist)],K)
    outres[i]<-mean(mn$y)
  }
  return(outres)
}

### sieve approximation with Bernstein polynomial ###
sieveB <- function(K,y,x){
  B <- rep(0,length(x))
  sieveB_dta <-as.data.table(cbind(y,1))
  z <- (x+2)/4 # normalizing x
  
  # set Bern. polynomial terms and solve
  for (i in 1:K){
    bk <- choose(K,i)*(z^i)*(1- z)^(K-i)
    sieveB_dta[,paste0("x",i):=bk]
  }
  Y<- as.matrix(y,ncol=1)
  X<- as.matrix(sieveB_dta[,y:=NULL])
  ols<- solve(t(X)%*%X, tol = 1e-20)%*%(t(X)%*%Y)
  out <- X%*%ols
  return(out)
}

### sieve approximation with linear splines ###
spline <- function(K,y,x){
  B <- rep(0,length(x))
  spline_dta <-as.data.table(cbind(y,1,x))
  qnum <- 1/(K+1) # quantile cutpoints
  
  # finding quantile cutpoints and estimating OLS within each
  spline_dta[,q:=as.numeric(cut(
    x,breaks=c(quantile(x, probs = seq(0,1,by=qnum))), labels=1:(1/qnum)))]
  spline_dta[is.na(q) & x <= -1.9, q :=1] # taking care of weird edge case
  for (i in 1:K){
    spline_dta[,paste0("x",i):= (q >= i)*(x - quantile(x,qnum*i))]
  }
  Y<- as.matrix(y,ncol=1)
  X<- as.matrix(spline_dta[,`:=`(q=NULL, y=NULL)])
  ols<- solve(t(X)%*%X, tol = 1e-20)%*%(t(X)%*%Y)
  out <- X%*%ols
  return(out)
}


################## Graphing Setup ################
# folder for saving figures
setwd(paste(filepath, "figures", sep="/"))
ty <- Y[,2]
tx <- X[,2]

# graphing tuning parameters
k_sieve<-10
bw_llr <- 0.2
bw_lcr <- 0.2

############## Monte Carlo Results ###############

# initializing input/output matrices
plotlist <- list()
methlist <- c("lcr", "llr", "sieve", "sieveB", "spline", "nn")
tunelist <- list(c(.02,.8), c(.1,.8) ,c(5,21),c(5,21),c(10,30),c(4,100))
argslist <- list('kern="unif"','kern="unif"'," ", " ", " ", 'match="maha"')
namelist <- list("Local Constant", "Local Linear", "Sieve", 
                 "Bern. Sieve", "Spline", "Nearest Neighbors")
tnamelist <- list("h", "h", "K", "K", "K", "K")
counter = 1

# looping over each nonparametric method
for (method in methlist){
  idx <- which(methlist == method)
  
  # looping over tuning parameters
  for (tuner in tunelist[[idx]]){
    
    # initializing data table for output
    assign(paste0("MC_",method),as.data.table(tx))
    
    # running functions over data (parallelized)
    ret<-foreach (d = 1:D) %dopar%{
      # extra arguments for specific methods
      if (method == "lcr" | method == "llr" | method =="nn"){
        rtn<-try(eval(parse(text=paste0(
          method,"(",tuner,", Y[,d], tx,",argslist[idx],")"))))
        if(class(rtn)=="try-error"){rtn<-NA} # skip singularity errors
      }
      if (method == "sieve" | method == "sieveB" | method=="spline") {
        rtn<-try(eval(parse(text=paste0(
          method,"(",tuner,", Y[,d], tx)"))))
        if(class(rtn)=="try-error"){rtn<-NA} # skip singularity errors
      }
      return(rtn)
    }
    
    # combining output with X values
    out <- do.call("cbind",ret)
    p <- cbind(get(paste0("MC_",method)),out)
    bGL <- as.data.table(melt(p, id=c("tx")))
    
    # calculating mean and sd of the estimates
    bGL[,sd:=sd(value),by=tx]
    t<-bGL[,lapply(.SD,mean),.SDcols=c("value","sd"), by=tx]
    names(t)<-c("x","mean", "sd")
    t[,`:=`(sdmin=mean-sd, sdmax=mean+sd)]
    
    # merging in average Y values across draws
    ytab <- as.data.table(cbind(tx,Y))
    ytab <- as.data.table(melt(ytab,id=c("tx")))
    s<-ytab[,mean(value), by=tx]
    setnames(s,"tx","x")
    t<-as.data.table(merge(t, s, by="x"))
    t<-t[order(x)]
    
    # plotting
    methname<-namelist[[idx]]
    tunename<-tnamelist[[idx]]
    plotlist[[counter]]<- ggplot(t, aes(x=x, y=mean)) + 
      geom_line(aes(x=x, y=V1), color="darkgoldenrod3")+ 
      geom_line() + 
      geom_ribbon(aes(ymin=sdmin, ymax=sdmax),
                  alpha=.65,fill="cornflowerblue", linetype=0) +
      theme_bw() + ggtitle(paste0(methname, ": ", tunename," = ",tuner)) + 
      xlab("X") + ylab("Y") +
      theme(plot.title = element_text(hjust = 0.5), 
            panel.grid.major = element_blank(), panel.grid.minor = 
              element_blank(), panel.background = element_blank(),
            axis.line = element_line(colour = "black"), 
            text = element_text(size = 12, family = myFont))
    
    # prep for next loop
    counter = counter + 1
    print(paste0("method ", method, " is done for tuning parameter ", 
                  tuner))
    
  }
}

# plotting graphs in two separate figures
graph1 <- grid.arrange(plotlist[[1]], plotlist[[2]],plotlist[[3]],
                       plotlist[[4]],plotlist[[11]],plotlist[[12]], ncol=2)
ggsave(graph1, width = 9, height = 11, file = file.path("tr_localnn.png"))

png("tr_sieves.png",width = 600, height = 800)
graph2 <- grid.arrange(plotlist[[5]], plotlist[[6]],plotlist[[7]],
                          plotlist[[8]],plotlist[[9]],plotlist[[10]],ncol=2)
ggsave(graph2, width = 9, height = 11, file = file.path("tr_sieves.png"))


# releasing cores for parallelization
stopImplicitCluster()


############## Plots by Kernel ###############
lcrout <-lcr(bw_lcr,ty,tx,"unif")
lcrout1 <-lcr(bw_lcr,ty,tx,"epan")
lcrout2 <-lcr(bw_lcr,ty,tx,"gauss")
lcrG <- as.data.table(cbind(tx,lcrout, lcrout1, lcrout2))
names(lcrG)<-c("X", "Uniform", "Epanechnikov", "Gaussian")
lcrGL <- as.data.table(melt(lcrG, id=c("X"))) #lcrGL = lcr output in long format
ytab <- as.data.table(cbind(tx,ty))
names(ytab) <- c("X","y")
lcrGL<-merge(lcrGL, ytab, by="X")

cc <- c("coral3", "cornflowerblue", "darkgoldenrod1")
ss <- c(.6,.6,.6)
lt <- "Legend"
png("tr_lcr_verify.png",width = 600, height = 450)
ggplot(data=lcrGL, aes(x=X, y=value, color=variable, size=variable)) +
  geom_point(aes(x=X, y=y), color="gray82") +
  geom_line() + theme_bw() + ylab("Y") + xlab("X") +
  scale_color_manual(lt,values=cc) + scale_size_manual(lt,values=ss) +
  theme(legend.title=element_text(size=13), legend.text=element_text(size=13)) +
  theme(plot.title = element_text(hjust = 0.5),
  text = element_text(size = 12, family = myFont)) +
  ggtitle("Local Constant")
dev.off()


### Local Linear Regression ###

# Plot by Kernel
llrout <-llr(bw_llr,ty,tx,"unif")
llrout1 <-llr(bw_llr,ty,tx,"epan")
llrout2 <-llr(bw_llr,ty,tx,"gauss")
llrG <- as.data.table(cbind(tx, llrout, llrout1, llrout2))
names(llrG)<-c("X", "Uniform", "Epanechnikov", "Gaussian")
llrGL <- melt(llrG, id=c("X")) #lcrGL = lcr output in long format
ytab <- as.data.table(cbind(tx,ty))
names(ytab) <- c("X","y")
llrGL<-merge(llrGL, ytab, by="X")

cc <- c("coral3", "cornflowerblue", "darkgoldenrod1")
ss <- c(.6,.6,.6)
lt <- "Legend"
png("tr_llr_verify.png",width = 600, height = 450)
ggplot(data=llrGL, aes(x=X, y=value, color=variable, size=variable)) +
  geom_point(aes(x=X, y=y), color="gray82") +
  geom_line() + theme_bw() + ylab("Y") + xlab("X") +
  scale_color_manual(lt,values=cc) + scale_size_manual(lt,values=ss) +
  theme(legend.title=element_text(size=13), legend.text=element_text(size=13))+
  theme(plot.title = element_text(hjust = 0.5),
  text = element_text(size = 12, family = myFont)) +
  ggtitle("Local Linear")
dev.off()





