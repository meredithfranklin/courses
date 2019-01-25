### R code from vignette source 'st_geos.Rnw'

# This code uses the gstat Irish wind data to illustrate the use of spatio-temporal statistics
library(spacetime) # for s-t analysis
library(gstat) # for s-t analysis (krigST)
library(mapdata) # for plotting map boundaries
library(maptools)
library(xts) # for time series

data(wind)
wind.loc$y = as.numeric(char2dms(as.character(wind.loc[["Latitude"]])))
wind.loc$x = as.numeric(char2dms(as.character(wind.loc[["Longitude"]])))
coordinates(wind.loc) = ~x+y
proj4string(wind.loc) = "+proj=longlat +datum=WGS84"

wind.loc

# Plot map of Irish wind data

plot(wind.loc, xlim = c(-11,-5.4), ylim = c(51,55.5), axes=T, col="red",cex.axis =.7)
text(coordinates(wind.loc), pos=1, label=wind.loc$Station, cex=.7)


# Date transformation
wind$time = ISOdate(wind$year+1900, wind$month, wind$day)
wind$jday = as.numeric(format(wind$time, '%j'))
stations = 4:15

# knots -> m/s
windsqrt = sqrt(0.5148 * as.matrix(wind[stations]))

# trend removal
Jday = 1:366
daymeans = sapply(split(windsqrt, wind$jday), mean)

meanwind = lowess(daymeans ~ Jday, f = 0.1)$y[wind$jday]
velocities = apply(windsqrt, 2, function(x) { x - meanwind })


# order locations to order of columns in wind, connect station names to location coordinates
wind.loc = wind.loc[match(names(wind[4:15]), wind.loc$Code),]
pts = coordinates(wind.loc[match(names(wind[4:15]), wind.loc$Code),])
rownames(pts) = wind.loc$Station
pts = SpatialPoints(pts)
# project to utm zone 29, to be able to do interpolation in
# proper Euclidian (projected) space:
proj4string(pts) = "+proj=longlat +datum=WGS84"
library(rgdal)
utm29 = CRS("+proj=utm +zone=29 +datum=WGS84")
pts = spTransform(pts, utm29)
# create space-time object
# constructed from space-wide dataset (each weather station is a column)
wind.data = stConstruct(velocities, space = list(values = 1:ncol(velocities)), time = wind$time, SpatialObj = pts)


# project the map boundary
m = map2SpatialLines(map("worldHires", xlim = c(-11,-5.4), ylim = c(51,55.5), plot=F))
proj4string(m) = "+proj=longlat +datum=WGS84"
m = spTransform(m, utm29)

# set up grid for predicting
grd = SpatialPixels(SpatialPoints(makegrid(m, n = 300)),proj4string = proj4string(m))
# select april 1961
wind.data = wind.data[, "1961-04"]

# select 10 prediction time points, evenly spread over this month
n = 10
tgrd = seq(min(index(wind.data)), max(index(wind.data)), length=n)
pred.grd = STF(grd, tgrd)

# S-T Variograms and Kriging
# separable covariance model, exponential with ranges 750 km and 1.5 day:
v = vgmST("separable", space = vgm(1, "Exp", 750000), time = vgm(1, "Exp", 1.5 * 3600 * 24), sill=0.6)

wind.ST = krigeST(values ~ 1, wind.data, pred.grd, v)
colnames(wind.ST@data) <- "sqrt_speed"

# Plot kriging results
layout = list(list("sp.lines", m, col='grey'),list("sp.points", pts, first=F, cex=.5))

print(stplot(wind.ST, col.regions=bpy.colors(), par.strip.text = list(cex=.5), sp.layout = layout))


# Plot time series at each station
library(lattice)
library(RColorBrewer)
b = brewer.pal(12,"Set3")
par.settings = list(superpose.symbol = list(col = b, fill = b), superpose.line = list(col = b),
        fontsize = list(text=9))

print(xyplot(values~time, groups=sp.ID, as.data.frame(wind.data), type='l', auto.key=list(space="right"),
        xlab = "1961", ylab = expression(sqrt(speed)), par.settings = par.settings))

# Example 3-D variogram (different dataset)
separableModel <- vgmST("separable", space=vgm(0.86, "Exp", 476, 0.14), time =vgm(   1, "Exp",   3, 0), sill=113)
data(vv)

if(require(lattice)) {

  plot(vv, separableModel, wireframe=TRUE, all=TRUE)

}

##### Using GAM models
library(mgcv)
library(tidyverse)
# spatio-temporal gam
# separable
# need data in long format

wind_long <- gather(wind, site, wind, RPT:MAL, factor_key=TRUE)
wind_long

# add coordinates to wind data
locs<-data.frame(coordinates(wind.loc)) %>% slice(rep(1:n(), each = dim(wind)))

wind_long<-data.frame(wind_long,locs)

sep.gam<-gam(wind~s(x,y,k=10,bs="ts")+s(month,k=8,bs="cr",fx=TRUE)+s(day,k=15,bs="cr",fx=TRUE),data=wind_long)
plot(sep.gam)

# example with tensor product can try with wind data
n <- 500
v <- runif(n)
w<-runif(n)
u<-runif(n)

test2<-function(u,v,w,sv=0.3,sw=0.4)  
{ ((pi**sv*sw)*(1.2*exp(-(v-0.2)^2/sv^2-(w-0.3)^2/sw^2)+
                  0.8*exp(-(v-0.7)^2/sv^2-(w-0.8)^2/sw^2)))*(u-0.5)^2*20
}

f <- test2(u,v,w)
y <- f + rnorm(n)*0.2

b <- gam(y~t2(v,w,u,k=c(25,5),d=c(2,1),bs=c("tp","cr"),full=TRUE),
         method="ML")
## more penalties now. numbers in labels like "r1" indicate which 
## basis function of a null space is involved in the term. 
pen.edf(b) 

# 3D plot
vis.gam(b)

