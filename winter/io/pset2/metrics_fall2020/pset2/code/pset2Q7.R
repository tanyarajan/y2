#####################################################
#               Metrics Problem Set 2
#                      Question 7
#                     Tanya Rajan
# Description: Reproducing Table 2 
#####################################################

# setup
rm(list=ls())
require(data.table)
require(ggplot2)
require(pracma)
require(knitr)
require(tidyr)
require(kableExtra)
require(foreign)
require(maxLik)
require(future.apply)
require(parallel)


# setting up parallelization
cor<-floor(detectCores(all.tests=FALSE)*.7)
if (is.na(cor) | cor==0){cor <- 1}
plan(multisession, workers = cor)

# set working directory (change to run)
filepath <- "/Users/tanyarajan/Documents/git/psets/metrics_fall2020/pset2" 
setwd(filepath) 

# read in data
datain <- as.data.table(read.csv("data/abadie.csv"))

# read in code file defining functions
source("code/pset2fns.R")


########### Column 1: OLS ########### 
# setting up regressors
y <- datain$nettfa
datain[,agemin := age - 25]
datain[,agemin2:=agemin^2]
x <- cbind(datain$p401k, datain$inc, datain$agemin, 
           datain$agemin2, datain$marr,datain$fsize)

# estimation and output
col1 <- reg(y*1000,x, se="HC1")
colnames <- c("constant", "p401k", "inc", "agemin", 
              "agemin2", "marr", "fsize")
c1out <- as.data.table(cbind(colnames, as.vector(col1$b), as.vector(col1$se)))
names(c1out) <- c("colnames", "beta1", "se1")


########### Column 2: First stage ########### 
# setting up regressors
xfs <- datain$p401k
zfs <- cbind(datain$inc, datain$agemin, datain$agemin2, 
             datain$marr,datain$fsize, datain$e401k)

# estimation and output
col2 <- reg(xfs,zfs, se="HC1")
colnames <- c("constant", "inc", "agemin", "agemin2", 
              "marr", "fsize", "e401k")
c2out <- as.data.table(cbind(colnames, as.vector(col2$b), as.vector(col2$se)))
names(c2out) <- c("colnames", "beta2", "se2")


########### Column 3: 2SLS ########### 
# setting up regressors
x2sls<- cbind(datain$p401k, datain$inc, datain$agemin, 
              datain$agemin2, datain$marr, datain$fsize)

# estimation and output
col3 <- iv(y*1000,x2sls,zfs, se="HC1")
colnames <- c("constant", "p401k", "inc", "agemin", 
              "agemin2", "marr", "fsize")
c3out <- as.data.table(cbind(colnames, as.vector(col3$b), as.vector(col3$se)))
names(c3out) <- c("colnames", "beta3", "se3")

  
# Anderson-Rubin confidence interval
testpar <- col3$b[2]
range <- 5000
mesh <- 300 # set low to allow the code to run faster
grid <- as.matrix(linspace(testpar-range, testpar+range, mesh), ncol=1)
ar.output <- future_apply(grid, 1, FUN=artest, betas=col3$b, bidx=2, y=y*1000, x=x2sls,z=zfs)
ar.ci <- c(grid[min(which(ar.output == 1))], grid[max(which(ar.output == 1))])


########### Column 4/5: LARF ########### 
# setting up regressors
xlarf<- cbind(datain$inc, datain$agemin, datain$agemin2, 
              datain$marr, datain$fsize)
marrlevels = unique(datain$marr)
agelevels = unique(datain$agemin)
c = 1
for (m in marrlevels){
  for (a in agelevels){
    datain[,paste0("dummy",c):= (marr == m & agemin==a)]
    c=c+1
  }
}
xdummies <- datain[,12:ncol(datain)] # dummies per marriage/age level

# sampling for bootstrapping standard errors
bsreps<-100
bootdta <- datain[,c("nettfa", "p401k", "e401k", 
                     "inc", "agemin", "agemin2", "marr", "fsize")]
bootdta[,nettfa:=nettfa*1000]
bootdta<-cbind(bootdta,xdummies)

bsamples<-bs(bootdta,bsreps)

### Attempt 1: nonparametric estimation and output (using K = 9 in sieve) ###
col4 <-larf(y=datain$nettfa*1000, d=datain$p401k, z=datain$e401k, 
            x=xlarf, xpred=cbind(datain$inc, xdummies), K=9, fs="nonparametric")
colnames <- c("constant", "p401k", "inc", "agemin", 
              "agemin2", "marr", "fsize")

# bootstrapping SEs
datalarf.np <- function(x){
  g <- as.data.table(x)
  return(larf(y=g$nettfa,d=g$p401k, z=g$e401k, 
              x=cbind(g$inc, g$agemin, g$agemin2, g$marr, g$fsize), 
              xpred=cbind(g$inc, g[,12:ncol(g)]), K=9, fs="nonparametric")$b)
}
est.np<-future_apply(bsamples,3,datalarf.np)
col4$se <- apply(est.np,1,sd)
c4out <- as.data.table(cbind(colnames, as.vector(col4$b), as.vector(col4$se)))
names(c4out) <- c("colnames", "beta4", "se4")


### Attempt 2: parametric estimation and output (using logit) ###
col5 <-larf(y=datain$nettfa*1000, d=datain$p401k, z=datain$e401k, 
            x=xlarf, xpred=cbind(datain$inc, xdummies), K=9, fs="parametric")
colnames <- c("constant", "p401k", "inc", "agemin", 
              "agemin2", "marr", "fsize")

# bootstrapping SEs
datalarf.p <- function(x){
  g <- as.data.table(x)
  return(larf(y=g$nettfa,d=g$p401k, z=g$e401k, 
              x=cbind(g$inc, g$agemin, g$agemin2, g$marr, g$fsize), 
              xpred=cbind(g$inc, g[,12:ncol(g)]), K=9, fs="parametric")$b)
}
est.p<-future_apply(bsamples,3,datalarf.p)
est.p <- est.p[,colSums(!is.na(est.p)) > 0]
col5$se <- apply(est.p,1,sd)
c5out <- as.data.table(cbind(colnames, as.vector(col5$b), as.vector(col5$se)))
names(c5out) <- c("colnames", "beta5", "se5")



########### Part c: Jackknife IV ##############

# main estimates
col6 <- iv(y*1000,x2sls,zfs, se="HC1", est="jack")
colnames <- c("constant", "p401k", "inc", "agemin", 
              "agemin2", "marr", "fsize")

# bootstrapping
bsreps <- 100

# drawing samples
bootdta <- datain[,c("nettfa", "p401k", "e401k", 
                     "inc", "agemin", "agemin2", "marr", "fsize")]
bootdta[,nettfa:=nettfa*1000]
bsamp.jack <- bs(bootdta, bsreps)

# function to parse data
datajack <- function(x){
  g <- as.data.table(x)
  return(iv(Y=g$nettfa, Z=cbind(g$inc, g$agemin, g$agemin2, 
                                g$marr,g$fsize, g$e401k), 
              X=cbind(g$p401k, g$inc, g$agemin, g$agemin2, g$marr, g$fsize), 
              est="jack", se="HC1")$b)
}


# running over samples
bsest.jack<-future_apply(bsamp.jack,3,datajack)
apply(as.matrix(bsest.jack), 1, sd)
c6out <- as.data.table(cbind(colnames, as.vector(col6$b), as.vector(col6$se)))
names(c6out) <- c("colnames", "beta6", "se6")


########### Outputting to a Table ##############
### Main Table ###
# combining all estimates
xfull <- c("constant", "p401k", "inc", "agemin", 
           "agemin2", "marr", "fsize", "e401k")
a1<-merge(c1out,c2out,"colnames", all.x=TRUE, all.y=TRUE)
a2<-merge(a1, c3out, by="colnames", all.x=TRUE, all.y=TRUE)
a3<-merge(a2, c4out, by="colnames", all.x=TRUE, all.y=TRUE)
ftab <- merge(a3, c5out, by="colnames",all.x=TRUE, all.y=TRUE)

# column order
ftab[colnames=="p401k",id := 1][colnames=="constant",id := 2][colnames=="inc", id:=3]
ftab[colnames=="agemin", id := 4][colnames=="agemin2", id:=5][colnames=="marr", id:=6]
ftab[colnames=="fsize", id:=7][colnames == "e401k", id:=8]
ftab<-as.data.table(apply(ftab[,-1], c(1,2), as.numeric))
ftab<-as.data.table(ftab[order(id)])

# format numbers
ftab<-as.data.table(apply(ftab, c(1,2), function(x) format(round(x,5), 
                                           digits=3, big.mark=",", nsmall=3)))

# create output table in correct format
outtable <- c(NA, NA, NA, NA, NA)
for (i in xfull){
  m <- which(xfull==i)
  outtable <- rbind(outtable, c(ftab$beta1[m], 
                                ftab$beta2[m], 
                                ftab$beta3[m], 
                                ftab$beta4[m],
                                ftab$beta5[m]))
  outtable <- rbind(outtable, c(as.character(paste0("(",ftab$se1[m],")")), 
                                as.character(paste0("(",ftab$se2[m],")")), 
                                as.character(paste0("(",ftab$se3[m],")")), 
                                as.character(paste0("(",ftab$se4[m],")")),
                                as.character(paste0("(",ftab$se5[m],")"))))
}
outtable<- outtable[-1,]
outtable[outtable=="NA"]<-" "
outtable[outtable=="(NA)"]<- " "

# row names
names <- c( "Participate in 401(k)", " ", "Constant", " ", "Family Income", " ",
           "Age", " ", "Age squared", " ", "Married", " ", "Family Size", " ",
           "Eligibility for 401(k)", " ")
outtable <- cbind(names, outtable)
colnames(outtable)<-c(" ", "OLS","First stage","Second Stage",
                      "LST (Nonpar)", "LST (Par)")

# write to tex
kable(outtable, format="latex", booktabs=TRUE) %>%
  add_header_above(c(" "=3, "2SLS" = 2, " " = 2)) %>%
  add_header_above(c(" "=3, "Endogenous treatment" = 4)) %>%
  write(file="tables/tr_table2.tex")

###### IV Table ######
xiv <- c("constant", "p401k", "inc", "agemin", 
           "agemin2", "marr", "fsize")
ftab <- merge(c3out,c6out,"colnames", all.x=TRUE, all.y=TRUE)
ftab[colnames=="p401k",id := 1][colnames=="constant",id := 2][colnames=="inc", id:=3]
ftab[colnames=="agemin", id := 4][colnames=="agemin2", id:=5][colnames=="marr", id:=6]
ftab[colnames=="fsize", id:=7][colnames == "e401k", id:=8]
ftab<-as.data.table(apply(ftab[,-1], c(1,2), as.numeric))
ftab<-as.data.table(ftab[order(id)])

# format numbers
ftab<-as.data.table(apply(ftab, c(1,2), function(x) format(round(x,5), 
                                                           digits=3, big.mark=",", nsmall=3)))

# create output table in correct format
outtable <- c(NA,NA)
for (i in xiv){
  m <- which(xiv==i)
  outtable <- rbind(outtable, c(ftab$beta3[m],ftab$beta6[m]))
  outtable <- rbind(outtable, c(as.character(paste0("(",ftab$se3[m],")")),
                                as.character(paste0("(",ftab$se6[m],")"))))
}
outtable<- outtable[-1,]
outtable[outtable=="NA"]<-" "
outtable[outtable=="(NA)"]<- " "

# row names
names <- c( "Participate in 401(k)", " ", "Constant", " ", "Family Income", " ",
            "Age", " ", "Age squared", " ", "Married", " ", "Family Size", " ")
outtable <- cbind(names, outtable)
colnames(outtable)<-c(" ", "TSLS","Jackknife TSLS")

# write to tex
kable(outtable, format="latex", booktabs=TRUE) %>%
  write(file="tables/tr_q7c.tex")

