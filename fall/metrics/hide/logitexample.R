####### Logit Function #######

require(maxLik)

###### Read in Data #########
# set working directory (change to run)
filepath <- "/Users/tanyarajan/Documents/git/psets/metrics_fall2020/pset3" 
setwd(filepath) 

# read in data
datain <- as.data.table(read.csv("data/angrist_evans_clean.csv"))
datain[,seq := seq(1:nrow(datain))]
datain <- datain[seq < 1000]

# variables of interest
Y <- as.matrix(datain$worked, ncol=1)
D <- as.matrix(datain$more2kids, ncol=1)
X <- as.matrix(datain[,list(age, ageat1st, agekid1, agekid2, boy1st, boy2nd, black, hispanic, otherrace)])


#### Logit Likelihood Function ###
fnL <- function(theta, yvec,xvec){
  z <- xvec%*%theta
  v1 <- log(exp(z)/(1 + exp(z)))
  v0 <- log(1/(1 + exp(z)))
  inner <- yvec*v1 + (1-yvec)*v0
  return(sum(inner))
}

#### Maximizing Likelihood ###
Xlog <- cbind(1, X)
guess_theta = rep(.1, ncol(Xlog))

# my function
mylogit <- maxLik(fnL, start=guess_theta, yvec = D, xvec= Xlog)

# compare to R
rlogit <- glm(D ~ X, family = "binomial")





