#################################################
# Example R code for PM599 spatial Lecture 2	  #
# Geostatistics									                #
#################################################

library(geoR)
library(maps)
library(rgdal)

pm.dat<-read.csv("pmdata.csv")

# Create geodata objects for original and logged PM2.5 concentrations
pm.dat.geo<-as.geodata(pm.dat,coords.col=c(1,2), data.col=4)
pm.dat.geo.log<-as.geodata(pm.dat,coords.col=c(1,2),data.col=3)

# Plotting the geodata object is useful for exploratory analysis
pdf("PMexploratory.pdf")
plot(pm.dat.geo)
dev.off()

pdf("PMexploratorylog.pdf")
plot(pm.dat.geo.log)
dev.off()

# Plotting as point data with different colour ramps
points(pm.dat.geo,cex.max=0.8,col=gray(seq(1, 0.1, l=100)), pt.divide="equal")

points(pm.dat.geo,cex.max=0.8,col= colorRampPalette(c("darkblue", "white", "darkred"))(20), pt.divide="equal")

points(pm.dat.geo,cex.max=0.8,col= rev(heat.colors(20)), pt.divide="equal", main="PM2.5 concentrations (ug/m3)", xlab="Longitude",ylab="Latitude")
# Add state boundaries
map("state",add=TRUE)

# subset geodata to look at a particular region
pdf("NY_PM25.pdf")
pm.dat.geo.ss<-subset(pm.dat.geo,lon> -75 & lat<41)
points(pm.dat.geo.ss,cex.max=1.5,col= rev(heat.colors(20)), pt.divide="equal", main="PM2.5 concentrations (ug/m3)", xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)
dev.off()

# projecting lat and lon to utm plane
# project() is from rgdal or proj4 library
newcoords<-project(as.matrix(cbind(pm.dat$lon, pm.dat$lat)), "+proj=utm +zone=18 +north")
pm.dat$x<-newcoords[,1]
pm.dat$y<-newcoords[,2]

# projecting with albers equal area
# mapproject() is from maptools library

xy<-mapproject(pm.dat$lon, pm.dat$lat, proj='albers', param=c(29.5,45.5), orientation=c(90,0,-80))
pm.dat$xx<-xy$x
pm.dat$yy<-xy$y

pm.dat.geo.utm<-as.geodata(pm.dat,coords.col=c(5,6), data.col=4)
plot(pm.dat.geo.utm)

pm.dat.geo.albers<-as.geodata(pm.dat,coords.col=c(5,6), data.col=4)
plot(pm.dat.geo.albers)

points(pm.dat.geo.albers,cex.max=0.8,col= rev(heat.colors(20)), pt.divide="equal", main="PM2.5 concentrations (ug/m3)", xlab="Easting",ylab="Northing")
# add state boundaries, needs to be projected
map("state",c("Massachusetts","Maine","Kentucky","Virginia","West Virginia","Maryland","Delaware","Vermont","New Hampshire","New York",
              "Ohio","Illinois","Michigan","Wisconsin","Indiana","Pennsylvania","New Jersey"),add=TRUE,projection="albers",param=c(29.5,45.5))

# empirical semivariograms: cloud, bin, boxplot bin
# use UTM so distances are in meters
vario1<-variog(pm.dat.geo.utm,option="cloud")
vario1$max.dist
pdf("VariogCloud.pdf")
plot(vario1,xlab="Distance (h), meters")
dev.off()

vario2<-variog(pm.dat.geo.utm,uvec=seq(0,vario1$max.dist,l=50),option="bin")
pdf("VariogBin.pdf")
plot(vario2,xlab="Distance (h), meters")
dev.off()


vario3<-variog(pm.dat.geo.utm,uvec=seq(0,vario1$max.dist,l=20),bin.cloud=T)
pdf("VariogBox.pdf")
plot(vario3,bin.cloud=T,xlab="Distance (h)")
dev.off()

# Robust estimator "modulus" is preferred
vario2m<-variog(pm.dat.geo.utm,uvec=seq(0,vario1$max.dist,l=20),option="bin",estimator.type="modulus")
pdf("VariogBinModulus.pdf")
plot(vario2,xlab="Distance (h), meters")
dev.off()

vario3m<-variog(pm.dat.geo.utm,uvec=seq(0,vario1$max.dist,l=20),bin.cloud=T,estimator.type="modulus")
pdf("VariogBoxModulus.pdf")
plot(vario3,bin.cloud=T,xlab="Distance (h), meters")
dev.off()

# Theoretical semivariograms by curve fitting (eyeball method)

# Exponential
pdf("VariogExpFit.pdf")
par(mfrow=c(1,2))
plot(vario2,xlab="Distance (h)")
curve(5+30*(1-exp(-x/1500000)),from = 0 , to=vario1$max.dist, ylab="Semivariance", xlab="Distance (h)",ylim=c(0,50),add=T)
dev.off()


# Spherical
pdf("VariogSphFit.pdf")
par(mfrow=c(1,2))
plot(vario2,xlab="Distance (h)")
curve(5+30*(1.5*(x/1500000)-0.5*(x/1500000)^3),from = 0 , to=vario1$max.dist, ylab="Semivariance", xlab="Distance (h)",add=T)
dev.off()

# Exploring anisotropy
vario.dir<-variog4(pm.dat.geo.utm,uvec=seq(0,vario1$max.dist,l=20),option="bin",estimator.type="modulus")
plot(vario.dir)
