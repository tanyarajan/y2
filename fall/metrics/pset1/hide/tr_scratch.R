# Old code for metrics pset 1


# Average lcr function
# I think this is wrong bc it returns a constant

lcr_avg <- function(h, y, x){
  lcr_dta <- as.data.table(cbind(y,x))
  outres <- rep(NA, length(x))
  for (i in 1:length(x) ){
    val = x[i]
    outres[i] <- lcr_dta[(val - x)^2 < h, mean(y)]
  }
  return(mean(outres)) # bias?
}

# test and plot
lcr(.001, Y[,1], X[,1])
test <- linspace(0,1, n=100)
ret <- sapply(test,FUN=function(x) lcr(x, Y[,1], X[,1]))
plot(test,ret)

ggplot() + xlim(0,5) + geom_function(fun = function(x) lcr(x, Y[,1], X[,1]))

bw <- npregbw(formula=Y[,1]~1, tol=.1, ftol=.1)
npreg(gradients=TRUE, bws=bw)


# Old LLR function  - gives funky results, think the derivation is wrong
llr <- function(h, y, x, kern){
  outres <- rep(NA, length(x))
  llr_dta <- as.data.table(cbind(1,y,x))
  for (i in 1:length(x)){
    val = c(1,x[i])
    iter <- cbind(llr_dta$V1, llr_dta$x)
    llr_dta[, dist:=apply(iter, 1, FUN=function(k) dist(rbind(val,k)))]
    llr_dta[, weight:=kernel(kern, (dist/h))]
    denom <- sum((llr_dta$x)^2 * llr_dta$weight)
    num <- sum(llr_dta$x * llr_dta$y * llr_dta$weight)
    outres[i]<- num/denom
  }
  return(outres)
}



# old monte carlo graphing code
system.time(){
  b<-apply(Y,c(1,2),FUN=lcr, h=.2,x=tx,kern="unif")
}

Y <- apply(X, c(1,2), function(x) sin(2*x)) + 
  apply(X, c(1,2), function(x) 2*exp(-16*(x^2))) + U


test<-apply(X, c(1,2), function(x) sin(2*x)) + 
  apply(X, c(1,2), function(x) 2*exp(-16*(x^2))) + U

brk<-50 # number of bins


# plotting means by bins 
bdta <- as.data.table(b)
bdta[,bins:=cut(bdta$V1, breaks=brk, labels=1:brk)]
bdta[,sd:=sd(V2), by=bins]
bdta[,binnum:=as.numeric(bins)]
t<-bdta[,lapply(.SD,mean),.SDcols=c("V2", "sd","binnum"), by=bins]
names(t)<-c("x","mean", "sd", "bins")
t[,`:=`(sdmin=mean-sd, sdmax=mean+sd)]
t[,var:=sd^2]

# Line plot with CI!
ggplot(t, aes(x=bins)) + 
  geom_line(aes(y=mean)) + 
  geom_ribbon(aes(ymin=sdmin, ymax=sdmax),alpha=.5,fill="cornflowerblue", linetype=0) +
  theme_bw()

# plotting variance separately
ggplot(t, aes(x=bins)) + 
  geom_line(aes(y=var)) 

# binning in graph sort of
bdta[,bins:=cut(bdta$V1, breaks=brk)]
bdta[order(V1)]
ggplot(bdta, aes(x=bins)) + stat_summary(aes(y=V2,group=bins), fun="mean",geom="point") + 
  stat_summary(aes(y=V2, group=bins), fun.data = mean_sdl, fun.args = list(mult=1), geom = "errorbar") +
  theme_bw()


# sd plot
m <- tapply(bdta$tries, bdta$bins, mean)
sd <- tapply(bdta$tries, bdta$bins, sd)
df <- data.frame(mean.y = m, sd = sd, bin = names(m))

ggplot(df, aes(x = bin, y = mean.y, 
               ymin = mean.y - 1.96*sd, 
               ymax = mean.y + 1.96*sd)) + 
  geom_errorbar() + geom_point(size = 3)


# eval/ parse function call stuff
assign(paste0("MC","test"),as.data.table(tx))

paste(c("t",methlist[[1]],10),collapse="")
fntest<-eval(parse(text=paste(c("lcr","(",tunelist[[1]][1],", Y[,2], tx,",argslist[1],")"),collapse="")))
fntest2<-eval(parse(text=paste(c("lcr","(",tunelist[[1]][1],", Y[,3], tx,",argslist[1],")"),collapse="")))

paste(c("lcr","(",tunelist[[1]][1],", ty, tx,",argslist[1],")"),collapse="")

tester <- eval(parse(text=fntest))

myoptions <- 'kern="unif"'
#foo <- eval(parse(text = paste("lcr(.2,ty,tx,", myoptions, ")")))
foo2 <- eval(parse(text = paste("lcr(.2,ty,tx,", myoptions, ")")))


# taking mean over function
max<-max(tx)
min<-min(tx)
mean <- (1/(max - min))*integrate.xy(tx,lcroutA)
sdev <- (1/N)*sqrt(sum((ty - lcrout)^2))


#--------------------------------------------------------------------

# random MC bullshit
system.time(
  ret<-foreach (d = 1:50) %dopar%{
    #tries <- sieve(40, Y[,1], X[,1])
    return(lcr(.02,Y[,d],tx,"unif"))
  }
)


cbind(MCtest,1)

cbind(get(paste0("MC","test")),1)

as.
out<-do.call("cbind", ret)
bGL <- as.data.table(melt(cbind(get(paste0("MC","test")),out),id=c("tx")))
bGL[,sd:=sd(value),by=tx]
t<-bGL[,lapply(.SD,mean),.SDcols=c("value","sd"), by=tx]
names(t)<-c("x","mean", "sd")
t[,`:=`(sdmin=mean-sd, sdmax=mean+sd)]

ytab <- as.data.table(cbind(tx,ty))
names(ytab) <- c("x","y")
t<-as.data.table(merge(t, ytab, by="x"))

t<-t[order(x)]
ggplot(t, aes(x=x, y=mean)) + 
  geom_point(aes(x=x, y=y), color="gray82")+ 
  geom_line() + 
  geom_ribbon(aes(ymin=sdmin, ymax=sdmax),alpha=.5,fill="cornflowerblue", linetype=0) +
  theme_bw()



### old monte carlo dummy loop ###
# tuning parameters
hlist <- c(0.1, 0.8)
Klist <- c(10, 20)

# looping over draws and tuning parameters 
# try to do apply to the draws? 
ktype<-"unif"
for (d in 1:D){ 
  for (t in 1:length(hlist)){
    lcr(hlist[t], Y[d], X[d], ktype)
    llr(hlist[t], Y[d], X[d], ktype)
  }
  for (t in 1:length(Klist)){
    #sieve(Y[d], X[d], Klist[t])
    #nn(Y[d], X[d], Klist[t])
    #sieveB(Y[d], X[d], Klist[t])
    #spline(Y[d], X[d], Klist[t])
  }
}

V2<-lcr(.02,Y[,1],X[,1],"unif")
l<-cbind(X[,1], V2)
b<-rbind(b,l)



