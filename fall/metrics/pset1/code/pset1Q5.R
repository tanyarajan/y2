#####################################################
#               Metrics Problem Set 1
#                      Question 5
#                     Tanya Rajan
# Description: This code cleans data 
#####################################################

# setup
rm(list=ls())
require(data.table)
require(ggplot2)
require(Hmisc)
require(foreign)
require(maxLik)
require(future.apply)

# set working directory (change to run)
filepath <- "/Users/tanyarajan/Documents/git/psets/metrics_fall2020/pset1" 
tables <- paste0(filepath,"/tables")
setwd(filepath) 

#################### Data Cleaning ###################
# load data
datain <- as.data.table(read.dta("data/Dataset_QJE_Replicate_with_Cities.dta"))
wdta <- datain[,c("pog1349", "c25prot", "c25juden", "pog20s", 
                "kreis", "kreis_nr", "c25pop", "judaica", "comm1349",
                "Latitude", "Longitude")]

# data cleaning
wdta[,`:=` (exist1349=(judaica == 1 | comm1349 == 1),
            logpop=log(c25pop), pProt=100*(c25prot/c25pop), 
            pJew=100*(c25juden/c25pop))]

# keeping only observations with pogroms in 1349
wdta<-wdta[exist1349==1]
wdta<-na.omit(wdta)

################ Method 1: Regression #################
# preparing variables
X <- as.matrix(wdta[,c("pog1349", "logpop",  "pJew", "pProt")])
Xols <- cbind(1,X)
Y <- as.matrix(wdta$pog20s, ncol=1)

# regression
ols <- solve(t(Xols)%*%Xols, tol=1e-20)%*%(t(Xols)%*%Y)
wdta[,yhat:=Xols%*%ols]

# clustered standard errors
klist <- unique(wdta$kreis_nr)
omega <- matrix(0, dim(Xols)[2], dim(Xols)[2])
for (g in klist){
  gdta <- wdta[kreis_nr == g,]
  xg <- cbind(1, as.matrix(gdta[,c("pog1349", "logpop",  "pJew", "pProt")]))
  epg <- gdta$pog20s - gdta$yhat
  omega <- omega + t(xg)%*%epg%*%t(epg)%*%xg
}
covmat <- solve(t(Xols)%*%Xols, tol=1e-20)%*%(omega)%*%solve(t(Xols)%*%Xols, tol=1e-20)
se <- sqrt(diag(covmat))

# adjusted R^2
tss <- sum((Y-mean(Y))^2)
ssr <- sum((Y - wdta$yhat)^2)
r2 <- 1 - ssr/tss
N <- length(Y)
adjR <- 1 - ((1-r2)*(N-1))/(N - ncol(X) - 1)

################ Method 2: NN Matching ################

# matching function
nnmatch <- function(K,y,x){
  
  # initializing
  obs<-nrow(x)
  matchid <- matrix(0, obs, obs)
  outres <- matrix(NA, obs,2)
  w<-var(x[,-1]) # calculate variance for mahalanobis
  
  # data frames
  nn_dta <- as.data.table(cbind(y,seq(1:obs),x))
  setnames(nn_dta,"V2","id")
  setnames(nn_dta,3,"D")
  d_1 <- nn_dta[D==1,] # treated
  d_0 <- nn_dta[D==0,] # untreated
  
  # looping through observations
    for (i in 1:obs){
      xvec<-as.matrix(x[i,-1],ncol=1)
      
      # matches from other treatment status
      oper<- eval(parse(text=paste0("d_",1-x[i,1])))
      match<-as.data.table(cbind(
        oper, apply(oper[,c(-1,-2,-3)], 1, # mahalanobis distance
                    function(k) sqrt(t(xvec - k)%*%solve(w)%*%(xvec - k)))))
      dist<-match[order(V2)]$V2[K]         # finding distance of Kth match
      mn<-match[V2<=as.double(dist),]      # finding matches within that dist.
      outres[i,1]<-mean(mn$V1)             # mean of y-values for match
                
      # for each id in opposite status matches, fill matchid matrix
      idlist<-as.vector(mn$id)  
      for (j in idlist){
        matchid[i,j]<-1
      }
      
      # matches from own treatment status (for variance)
      own <- eval(parse(text=paste0("d_",x[i,1])))
      vmatch<-as.data.table(cbind(
        own, apply(own[,c(-1,-2,-3)], 1,   # mahalanobis distance
                    function(k) sqrt(t(xvec - k)%*%solve(w)%*%(xvec - k)))))
      vdist<-vmatch[order(V2)]$V2[K+1]     # finding distance of Kth match
      vmn<-vmatch[V2<=as.double(vdist),]   # finding matches within that dist.
      vmn[,meanval:=mean(V1)]              # mean of y-values for match
      vmn[,err:=(V1-meanval)^2]            # deviations
      
      # estimator for heteroskedastic sigma
      outres[i,2]<-sum((1/(length(vmn$err)-1))*vmn$err)
    }

  # calculating K: weights based on how often observations were matched
  numtrt <- nrow(d_1)
  trt <- x[,1]
  tk<-c()
  denom <- rowSums(matchid)
  for (i in 1:length(y)){
    tk[i]<- sum(matchid[,i]*(1/denom))
  }

  # variance estimator
  sigma<-outres[,2]
  v<-(1/(numtrt^2))*sum(((trt-(1-trt)*tk)^2)*sigma)
  
  return(cbind(outres,v))
}

# Panel B specification
t<-nnmatch(4,Y,X)

# constructing counterfactuals to calculate ATT
wdta[,`:=`(cf=t[,1], card=t[,2])]
wdta[pog1349==0,`:=`(yhat0=as.double(pog20s),yhat1=cf)]
wdta[pog1349==1,`:=`(yhat0=cf,yhat1=pog20s)]
wdta[pog1349==1,att:=pog20s - yhat0]
att <- wdta[pog1349==1]

# ATT and variance from Method 2
ATT2<-mean(att$att,na.rm=TRUE)
seATT2<-sqrt(t[1,3])

############### Method 3: Geo Matching ################
# preparing geographic covariates
Xgeo <- as.matrix(wdta[,c("pog1349", "Latitude","Longitude")])

# Panel C specification
g<-nnmatch(4,Y,Xgeo)

# constructing counterfactuals to calculate ATT
wdta[,`:=`(g_cf=g[,1], g_card=g[,2])]
wdta[pog1349==0,`:=`(g_yhat0=as.double(pog20s),g_yhat1=g_cf)]
wdta[pog1349==1,`:=`(g_yhat0=g_cf,g_yhat1=pog20s)]
wdta[pog1349==1,g_att:=pog20s - g_yhat0]
g_att <- wdta[pog1349==1]

# ATT and variance from Method 3
ATT3<-mean(g_att$g_att,na.rm=TRUE)
seATT3<-sqrt(g[1,3])



############### Method 4: PScore Match ################

# Logit Likelihood Function
fnL <- function(theta, yvec,xvec){
  z <- xvec%*%theta
  v1 <- log(exp(z)/(1 + exp(z)))
  v0 <- log(1/(1 + exp(z)))
  inner <- yvec*v1 + (1-yvec)*v0
  return(sum(inner))
}

# Propensity Score matching function
psmatch <- function(data){
  # Setting things up
  yvec <- data[,1]
  dvec <- data[,2]
  xvec <- data[,c(-1,-2)]
  
  # Logit
  Xlog <- cbind(1,xvec)
  thetat <- rep(.1,ncol(Xlog)) # starting guess
  logit<-maxLik(fnL, start=thetat, yvec=X[,1],xvec=Xlog)
  pscore<- exp(Xlog%*%logit$estimate)/(1 + exp(Xlog%*%logit$estimate))
  
  # Matching on pscores
  psmat <- cbind(dvec,pscore)
  psm<-nnmatch(4,Y,psmat)
  
  # constructing counterfactuals to calculate ATT
  psd <- as.data.table(cbind(yvec,dvec))
  names(psd)<-c("yvec","dvec")
  psd[,`:=`(cf=psm[,1])]
  psd[dvec==0,`:=`(yhat0=as.double(yvec), yhat1=cf)]
  psd[dvec==1,`:=`(yhat0=cf, yhat1=yvec)]
  psd[dvec==1,att:=yvec - yhat0]
  psd[,ate:=yhat1 - yhat0]
  psm_att <- psd[dvec==1]
  
  # Propensity Score ATT
  ps_ATT<-mean(psm_att$att,na.rm=TRUE)
  
  # Propensity Score ATE
  ps_ATE<-mean(psd$ate,na.rm=TRUE)
  
  return(c(ps_ATT, ps_ATE))
}

# Panel D coefficients
ps_est<-psmatch(cbind(Y,X[,1], X[,-1]))
ATT4<- ps_est[1]
ATE4<- ps_est[2]

  
#### Bootstrapping Standard errors ####
# setting up bootstrap data
set.seed(1234)
bsreps <- 1000
sampling <- wdta[,c("pog20s","pog1349","logpop","pJew","pProt")]
sampling[,id:=seq(1:length(Y))]

# function to draw bootstrap samples
bs <- function(data, reps){
  outlet <- array(dim=c(nrow(data), ncol(data), reps))
  for (r in 1:reps){
    bstest <- sample.int(nrow(data), nrow(data))
    bsout <- matrix(NA,1,ncol(data))
    for (l in bstest){
      bsout<-rbind(bsout, data[l,],use.names=FALSE)
    }
    bsout<-bsout[-1,]
    outlet[,,r]<-as.matrix(bsout)
  }
  return(outlet)
}

# drawing bootstrap samples
bsamples<-bs(sampling,bsreps)


# Panel D standard errors
plan(multisession)
est<-future_apply(bsamples,3,psmatch)
seATT4 <- sd(est[1,])
seATE4 <- sd(est[2,])

# writing to tex tables
outtable <- matrix(" ",nrow=23,ncol=1)
counter = 2
for (c in 1:ncol(X)){
  outtable[c*2 ]<-as.character(round(ols[c+1],5))
  outtable[c*2 + 1]<-as.character(paste0("(",round(se[c+1],5),")"))
  counter = counter + 2
}
addl <- c(N, adjR)
for (j in addl){
  outtable[counter]<- as.character(round(j,3))
  counter = counter + 1
}
counter = counter + 1
telist <- c("ATT2", "ATT3", "ATT4", "ATE4")
for (j in telist){
  te <- eval(parse(text=j))
  sej <- eval(parse(text=paste0("se",j)))
  outtable[counter] <- as.character(round(te,5))
  outtable[counter+1] <- as.character(paste0("(",round(sej,5),")"))
  counter = counter + 3
}
names <- c("Panel A", "Pogrom 1349", " ", "ln(Pop)", " ", 
           "%Jewish", " ", "%Protestant", " ", "N", "Adj. R2", 
           "Panel B", "Match ATT", " ", 
           "Panel C", "Geo Match ATT", " ", 
           "Panel D", "PS ATT"," ", "PS ATE", " ", " ")
outtable <- cbind(names, outtable)
                             
# saving table
write(kable(outtable,format="latex",booktabs=TRUE, 
            col.names=c(" ", "Estimates")),
      file=file.path(tables,"tr_replication.tex"))




