## Asymptotics
## PM522b Code for Slides 6

library(boot)

# Visualization of CLT
chi.mean<-NULL
for(i in 1:1000){
x<-rchisq(5,df=6)	
mx<-mean(x)
chi.mean<-c(chi.mean,mx)
}
plot(density(chi.mean),main="Random Sample Chi-Sq Distribution n=5",xlab="Observed Sample Mean",xlim=c(0,14))

chi.mean.large<-NULL
for(i in 1:1000){
x<-rchisq(25,df=6)	
mx<-mean(x)
chi.mean.large<-c(chi.mean.large,mx)
}

plot(density(chi.mean.large),main="Random Sample Chi-Sq Distribution n=25",xlab="Observed Sample Mean",xlim=c(0,14))


chi.mean.larger<-NULL
for(i in 1:1000){
x<-rchisq(150,df=6)	
mx<-mean(x)
chi.mean.larger<-c(chi.mean.larger,mx)
}

plot(density(chi.mean.larger),main="Random Sample Chi-Sq Distribution n=150",xlab="Observed Sample Mean",xlim=c(0,14))

# Simulation of mean of random standard normals with increasing sample size from 0 to 1000
x<-cumsum(rnorm(1000,mean=0,sd=1))/(1:1000)
plot(1:1000,x,type="l",xlab="Iteration",ylab="Average")
abline(h=0,lty=2)


###################################
# Asymptotic Likelihood Ratio Test
n = 50
p0 = 0.6
nSim = 10000
xPlotMax = 8
calcLRT <- function(xSum, n, p0) {
  pHat <- xSum/n
  return(-2*( xSum*(log(p0)-log(1-p0)-log(pHat)+log(1-pHat)) + n*(log(1-p0)-log(1-pHat))))
}
simLRT <- sapply(rbinom(nSim, n, p0), calcLRT, n=n, p0=p0)
pdf("LRTplot.pdf")
truehist(simLRT, xlim=c(0,xPlotMax), xlab="Likelihood Ratio Test Statistic Values",
         main="Simulated Statistics and Their Limiting Chi-squre Density",nbins=15)
xEvalPoints <- seq(from=0, to=xPlotMax, length.out=100)
lines(xEvalPoints, sapply(xEvalPoints, dchisq, df=1))
dev.off()


# Bootstrap
set.seed(1234)
mean.boot <- function(x,index){
 mean(x[index])
}
var.boot <- function(x,index){
  var(x[index])
}

# Non-parametric bootstrap
sample<-c(5,14,10,7)
boot.mean<-boot(sample, mean.boot, 1000)
mean(boot.mean$t)
var(boot.mean$t)

(1/4^4)*sum((boot.mean$t-mean(boot.mean$t))^2)

boot.var<-boot(sample,var.boot,1000)
mean(boot.var$t)
var(boot.var$t)

# Parametric bootstrap: sample is from a distribution

parametric.sample<-rnorm(2.71,4.82,n=9)
boot.mean.par<-boot(parametric.sample, mean.boot, 1000)


                                 