library(spatstat)
library(splancs)

#data bramblecanes,finpines,and spruces

# tree data
data(bramblecanes)
data(finpines)
data(spruces)


#creates point pattern data format from the x and y coordinates
finpines.pts<-as.points(finpines$x,finpines$y)
bramblecanes.pts <-as.points(bramblecanes$x,bramblecanes$y)
spruces.pts <-as.points(spruces$x,spruces$y)


#divide into quadrants and count number of points in the quadrant
Q <- quadratcount(finpines, nx = 4, ny = 3)
Q

plot(finpines.pts,pch=19, main="Finpines", xlab="Easting",ylab="Northing")
plot(Q, add = TRUE, cex = 2)


# Compare Finpines to CSR and regular process
sp_point <- matrix(NA, nrow=length(finpines$x),ncol=2)
sp_point[,1] <- finpines$x
sp_point[,2] <- finpines$y
colnames(sp_point) <- c("x","y")
plot(x=sp_point[,1],y=sp_point[,2],main="Finpines Data", xlab="Easting",ylab="Northing",cex=.5)

# Random points
u.x <- runif(n=nrow(sp_point), min=bbox(sp_point)[1,1], max=bbox(sp_point)[1,2])
u.y <- runif(n=nrow(sp_point), min=bbox(sp_point)[2,1], max=bbox(sp_point)[2,2])

# Regular points

r.x <- seq(from=min(sp_point[,1]),to=max(sp_point[,1]),length=sqrt(nrow(sp_point)))
r.y <- seq(from=min(sp_point[,2]),to=max(sp_point[,2]),length=sqrt(nrow(sp_point)))
r.x <- jitter(rep(r.x,length(r.x)),.001)
r.y <- jitter(rep(r.y,each=length(r.y)),.001)


par(mfrow=c(1,3),mar=c(4,4,1.5,0.5))
plot(x=sp_point[,1],y=sp_point[,2],main="Finpines Data", xlab="Easting",ylab="Northing",cex=.5,pch=19)
plot(x=u.x,y=u.y,main="Random Points", xlab="Easting",ylab="Northing",cex=.5,pch=19)
plot(x=r.x,y=r.y,main="Regular Points", xlab="Easting",ylab="Northing",cex=.5,pch=19)



# Complete spatial randomness

# CSR from random uniform point process
# Generates N uniform random points

par(mfrow=c(2,2),mar=c(0.4,0.4,1,0.4))
pp.100.unif<-runifpoint(100,win=owin(c(0,1),c(0,1)))
plot(pp.100.unif, main="")


# CSR from homogeneous poisson process

# CSR with intensity 100 in the unit square
pp.100.unit <- rpoispp(100,win=owin(c(0,1),c(0,1)))

# CSR with intensity 1 in a 10 x 10 square
pp.1.10 <- rpoispp(1, win=owin(c(0,10),c(0,10)))

par(mfrow=c(1,2))
plot(pp.100.unit,main="")
title('intensity = 100, unit square')
plot(pp.1.10,main="")
title('intensity = 1, 10 x 10 square')

# CSR with intensity 200 in the unit square
par(mfrow=c(1,1))
pp.200.unit <- rpoispp(200)

plot(pp.200.unit)

# Testing for CSR with Ripley's K
# create polygon (square) area for finpines data
minfx <- min(finpines$x)
maxfx <- max(finpines$x)
minfy <- min(finpines$y)
maxfy <- max(finpines$y)

polyfin <- as.points(c(minfx,maxfx,maxfx,minfx),c(minfy,minfy,maxfy,maxfy))
UL.khat.fin <- Kenv.csr(length(finpines$x), polyfin, nsim=99, seq(0.05,3.5,0.15))
khat.fin<-khat(finpines.pts, polyfin, seq(0.05,3.5,0.15))
# plot of Khat-pi t^2 from data
plot(seq(0.05,3.5,0.15), khat.fin-pi*seq(0.05,3.5,0.15)^2, type="l", xlab="Distance", ylab="Estimated K-pi h^2")

# plot upper bound
lines(seq(0.05,3.5,0.15), UL.khat.fin$upper-pi*seq(0.05,3.5,0.15)^2, lty=2)

# plot lower bound
lines(seq(0.05,3.5,0.15), UL.khat.fin$lower-pi*seq(0.05,3.5,0.15)^2, lty=2)

# Other plots (Variance stabilizing)
# L functions
l<-function(k,h)
{sqrt(k/pi)-h}

# plot Lhat from data
plot(seq(0.05,3.5,0.15), l(khat.fin, seq(0.05,3.5,0.15)), type="l", xlab="Distance", ylab="Estimated L",ylim=c(-.1,.2))

# plot upper bound of Lhat
lines(seq(0.05,3.5,0.15), l(UL.khat.fin$upper,seq(0.05,3.5,0.15)), lty=2)

# plot lower bound of Lhat
lines(seq(0.05,3.5,0.15), l(UL.khat.fin$lower,seq(0.05,3.5,0.15)), lty=2)

# different R function for Ripley's K
K<-Kest(finpines)
plot(K)
plot(K, cbind(r, sqrt(iso/pi)) ~ r)
plot(K, cbind(trans,iso,border) - theo ~ r)

plot(envelope(finpines,Kest))

# distance-based functions testing for CSR
Gtest<-Gest(finpines)
plot(envelope(finpines,Gest))

Ftest<-Fest(finpines)
plot(envelope(finpines,Fest))

Jtest<-Jest(finpines)
plot(envelope(finpines,Jest))

# Inhomogenous Poisson Process with different intensity functions

lamb1 <- function(x,y) {100 * exp(10*x-5*y)}
Ipp1 <- rpoispp(lamb1, 100)
Ipp1.p<- as.points(Ipp1$x,Ipp1$y)


plot(Ipp1.p,main="IPP, intensity(x,y)=100*exp(10x-5y)",pch=20,xlab='x',ylab='y')


lamb2 <- function(x,y) {100 * exp(-10*x+5*y)}
Ipp2 <- rpoispp(lamb2, 100)
Ipp2.p<- as.points(Ipp2$x,Ipp2$y)


plot(Ipp2.p,main="IPP, intensity(x,y)=100*exp(-10x+5y)",pch=20,xlab='x',ylab='y');


# can save a step and write function within rpoispp
Ipp3 <- rpoispp(function(x,y) {100 * exp(-3*x)}, 100)
plot(Ipp3,pch=20,xlab='x',ylab='y')


# Clutered point process
# Simulations of (2) clustered point processes
# Generate a square polygon

spoly<-as.points(c(0,1,1,0),c(0,0,1,1)) #polygon

#intensity of parents is 15, average number of offspring is 7, and variance of distance traveled by child (away from parent) in 0.0025
pcp1<-pcp.sim(15,7,0.00025,spoly)
plot(pcp1,main="PCP, (P,O,Spread)=(25,7,.00025)",pch=20,xlab='x',ylab='y')


#intensity of parents is 25, average number of offspring is 4, and variance of distance traveled by child in 0.005
pcp2 <-pcp.sim(25,4,0.005,spoly)
plot(pcp2,main="PCP, (P,O,Spread)=(25,4,.005)",pch=20,xlab='x',ylab='y')


# Estimating the parameters of a cluster process
# using the pcp2 data from above

# Using K-function for Clustered Process
# Simulate homogeneous pcp (same as above)
pcp2 <-pcp.sim(25,4,0.005,spoly)
# need to define polygonal area
spoly<-as.points(c(0,1,1,0),c(0,0,1,1))
# estimate
pcp.fit<-pcp(pcp2,spoly,h0=.3,n.int=40)

#simulate from the theoretical process with parameter estimates from above

plot(pcp2)
m.h <- npts(pcp2)/(areapl(spoly)*pcp.fit$par[2])
sims.h <- pcp.sim(pcp.fit$par[2], m.h, pcp.fit$par[1], spoly)
pointmap(as.points(sims.h), add=TRUE, col="red")

# L-function
l<-function(k,s)
{sqrt(k/pi)-s}

# K function for PCP
r <- seq(0.01,1,by=.2)
K.env <- Kenv.pcp(pcp.fit$par[2], m.h, pcp.fit$par[1], spoly, nsim=20, r=r)
L.env<- lapply(K.env,FUN=l,r)

limits <- range(unlist(L.env))
plot(r, sqrt(khat(as.points(pcp2),spoly,r)/pi)-r, ylim=limits, main="L function with simulation envelopes and average", type="l", xlab="distance", ylab="")
lines(r, L.env$lower, lty=5)
lines(r, L.env$upper, lty=5)
lines(r, L.env$ave, lty=6)
abline(h=0)

# Application to a real data set (Cardiff)
data(cardiff)
polymap(cardiff$poly)
pointmap(as.points(cardiff), add=TRUE)
title("Locations of homes of 168 juvenile offenders")
pcp.fit.cardiff <- pcp(as.points(cardiff), cardiff$poly, h0=20, n.int=30)
pcp.fit.cardiff

m <- npts(as.points(cardiff))/(areapl(cardiff$poly)*pcp.fit.cardiff$par[2])
r <- seq(2,30,by=2)
K.env <- Kenv.pcp(pcp.fit.cardiff$par[2], m, pcp.fit.cardiff$par[1], cardiff$poly,
                  nsim=20, r=r)
L.env<- lapply(K.env,FUN=l,r)

limits <- range(unlist(L.env))
plot(r, sqrt(khat(as.points(cardiff),cardiff$poly,r)/pi)-r, ylim=limits,
     main="L function with simulation envelopes and average", type="l",
     xlab="distance", ylab="")
lines(r, L.env$lower, lty=5)
lines(r, L.env$upper, lty=5)
lines(r, L.env$ave, lty=6)
abline(h=0)


