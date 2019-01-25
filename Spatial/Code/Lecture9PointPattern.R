
library(splancs)
library(spatstat)
library(dbscan)

#####################################
# Kernel density plots and tesslations
# A “tessellation” is a division of space 
# into non-overlapping regions. 
# Dirichlet(X) computes the Dirichlet
# tessellation or Voronoi tessellation of
# the point pattern X. 
######################################

# Clutered point process
spoly<-as.points(c(0,1,1,0),c(0,0,1,1)) #polygon
# intensity of parents is 25, average number of offspring is 4, and variance of distance traveled by child in 0.0025
pcp1<-pcp.sim(25,4,0.00025,spoly)
plot(pcp1)
# intensity of parents is 25, average number of offspring is 4, and variance of distance traveled by child in 0.005
pcp2 <-pcp.sim(25,4,0.005,spoly)
plot(pcp2)

# density and tesselation require as.ppp
pcp2.pp<-as.ppp(pcp2,c(0,1,0,1))
plot(density.ppp(pcp2.pp,0.1), main="Gaussian Kernel Density of PCP")


points(pcp2.pp, pch=19, cex=0.5)

plot(dirichlet(pcp2.pp), main="Dirichlet Tesselation of PCP")

# kernel density on data - use Spanish fires
plot(density.ppp(clmfires,10))
points(clmfires,pch=19,cex=0.1)


# kernel density and tesselation of leukemia data
plot(density.ppp(humberside)) 
points(humberside)
plot(dirichlet(humberside))



#####################################
# Fitting Point Processes
######################################

# Homogenous Poisson Process or CSR

# CSR with intensity 100 in the unit square
pp.100.unit <- rpoispp(100)

# Inhomogenous Poisson Process
# Simulate process using intensity lambda
lambda <- function(x,y) {exp(10+10*x-5*y)}
ipp1 <- rpoispp(lambda, 100)

# Fit inhomogeneous poisson process with ppm() function
ipp.fit1<-ppm(ipp1,~x+y)
summary(ipp.fit1)
plot(ipp.fit1) # this plots the density

# Simulate another IPP
ipp2 <- rpoispp(function(x,y) {100 * exp(-3*x)}, 100)

# Fit 
ipp.fit2<-ppm(ipp2,~x)
summary(ipp.fit2)
plot(ipp.fit2)
AIC(ipp.fit2)
     
#Clutered point process
spoly<-as.points(c(0,1,1,0),c(0,0,1,1)) #polygon

#intensity of parents is 25, average number of offspring is 4, and variance of distance traveled by child in 0.0025
pcp1<-pcp.sim(25,4,0.00025,spoly)
plot(pcp1)
#intensity of parents is 25, average number of offspring is 4, and variance of distance traveled by child in 0.005
pcp2 <-pcp.sim(25,4,0.005,spoly)
plot(pcp2)

# must put in ppp object (and specify window) for kppm to work
pcp2.pp<-as.ppp(pcp2,c(0,1,0,1))

# Fit Poisson Cluster Process to simulated data by minimum contrast
pcp.mod1<-kppm(pcp2.pp,method="mincon")
summary(pcp.mod1)

# Fit Log Gaussian Poisson cluster process to simulated data
lgcp.mod1<-kppm(pcp2.pp,clusters="LGCP",method="mincon")
summary(lgcp.mod1)
plot(lgcp.mod1) # is this on residuals?
plot(envelope(lgcp.mod1))

# Fit log Gaussian Poisson cluster process with matclust.estK function uses deviance
# use actual data (redwood)
lgcp.mod2 <- lgcp.estK(redwood, c(var=1, scale=0.1))
lgcp.mod2
plot(lgcp.mod2)

lgcp.redw<-kppm(redwood,clusters="LGCP",method="mincon")
summary(lgcp.redw)
plot(lgcp.redw)
plot(envelope(lgcp.redw))

# Fit Matern cluster process to simulated data
mat.mod1<-kppm(pcp2.pp,clusters="MatClust",method="mincon")
summary(mat.mod1)
plot(mat.mod1)
plot(envelope(mat.mod1))

# Fit Matern cluster process to redwood data
mat.mod2<-kppm(redwood,clusters="MatClust",method="mincon")
plot.kppm(mat.mod2)
summary(mat.mod2)
plot(envelope(mat.mod2))

# simulate data based on this fit
sim.matmod<-simulate(mat.mod2, 10)
# see what you get in each simulation
names(sim.matmod[[1]])
plot(sim.matmod)

# Fit Matern cluster process with inhomogenous intensity to redwood data
mat.mod3<-kppm(redwood ~x+y,clusters="MatClust",method="mincon")
summary(mat.mod3)

# Fit Matern cluster process with matclust.estK function uses deviance
mat.mod4 <- matclust.estK(redwood, c(kappa=10, scale=0.1))
mat.mod4
plot(mat.mod4)

# Inhibited Poisson Process (Regular pattern)
# Strauss models r is interaction distance
# Use biological cells point pattern data (42 cells on unit square)
data(cells)
plot(cells,pch=20)

# Homogenous Strauss process
fit1<-ppm(cells, ~1, Strauss(r = 0.1))
fit1
# Inhomogenous Strauss process
fit2<-ppm(cells, ~x+y, Strauss(r = 0.1))
fit2

# Inhomogenous Strauss process with more complicated log intensity
fit3<-ppm(cells, ~polynom(x, y, 2), Strauss(r = 0.1))
par(mfrow = c(1, 2))
plot(fit3, how = "image", ngrid = 256)

# Fit homogeneous Strauss process to Swedish pines data
data(swedishpines)
plot(swedishpines,pch=20)
fit4 <- ppm(swedishpines, ~1, Strauss(r = 7))
# Simulate data based on fit
sim.points <- rmh(fit4)
plot(sim.points, main = "Simulation from fitted Strauss model",pch=20)

# ppm fit diagnostics
diagnose.ppm(fit4)
qqplot.ppm(fit4, nsim = 39)

#####################################
# Density-based Clustering with DBSCAN 
# Hierarchical DBSCAN (HDBSCAN*)
# Using satellite data from VIIRS 
# Identifying clusters of flaring around 
# fracking sites
######################################

viirs <- read.csv('/Users/mf/Dropbox (University of Southern California)/Courses/PM569/data/viirs_2012.csv')
viirs.pts <- viirs[,5:6]

library(proj4)
proj4string <- "+proj=utm +zone=14 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs "

# Transformed data back to lat lon for leaflet plotting only
latlon <- project(viirs.pts, proj4string, inverse=TRUE)
viirs.latlon <- data.frame(lat=latlon$y, lon=latlon$x)

library(leaflet)
leaflet(viirs.latlon) %>% addProviderTiles("OpenStreetMap")  %>%
  addCircleMarkers(opacity = 1, radius=0.4, color = "blue") %>%
 addLegend("bottomright",  labels = "VIIRS Nightfire",colors="blue",opacity=1)

# subset to smaller region
viirs.pts<-viirs.pts[viirs.pts$x<500000,]

#dbscan
db.clust = dbscan(viirs.pts/1000, minPts = 5, eps=4) # what is epsilon (units?)
hullplot(viirs.pts/1000,db.clust)


hdb.clust = hdbscan(viirs.pts/1000, minPts = 3)
hullplot(viirs.pts/1000,hdb.clust)

# cluster membership assignment for each
hdb.clust$cluster
# cluster membership probability of each
hdb.clust$membership_prob
# outlier score for each point, *not* complement of membership scor
hdb.clust$outlier_scores
# cluster score for each *cluster*
hdb.clust$cluster_scores

# plot hierarchical tree
plot(hdb.clust)
plot(hdb.clust, show_flat = TRUE)

