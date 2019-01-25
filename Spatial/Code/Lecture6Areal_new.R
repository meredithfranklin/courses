library(spdep) # for areal data analysis
library(maptools) # for reading arcgis shapefiles
library(sf) #for reading arcgis shapefiles (new)
library(leaflet)
library(tidyverse)
library(pgirmess)


# Shapefile contains all of the data and polygon information
# the NC shapefile is also part of the sf library
# You can find it here (MacOS) /Library/Frameworks/R.framework/Versions/3.5/Resources/library/sf/shape/nc.shp

nc<-st_read(system.file("shape/nc.shp", package = "sf"), stringsAsFactors = FALSE)
plot(nc[10])


# Plot rate on leaflet
pal <- colorNumeric(palette = "YlGnBu",domain = nc$SID74/nc$BIR74)
nc<-st_transform(nc,4326) # need to transform so that shapefile has WGS
leaflet(nc) %>% addProviderTiles("OpenStreetMap")  %>%
  addPolygons(weight = 1, color = "#444444", opacity = 1,
              fillColor = ~pal(SID74/BIR74), fillOpacity = 0.9, smoothFactor = 0.5)%>%
  addLegend("bottomright", pal = pal, values = ~SID74/BIR74, title = "SIDS74 rate", opacity = 1)

# convert sf to spatial to use spdep
nc_sp <- sf:::as_Spatial(nc$geom)

plot(nc_sp)
#text(nc$FIPSNO)
# k nearest neighbours
# knearneigh() creates matrix with index for the regions belonging to knn
# knn2nb() creates neighbourhood list
sids_kn1<-knn2nb(knearneigh(coordinates(nc_sp), k=1, RANN=FALSE))
sids_kn2<-knn2nb(knearneigh(coordinates(nc_sp), k=2, RANN=FALSE))
sids_kn4<-knn2nb(knearneigh(coordinates(nc_sp), k=4, RANN=FALSE))

# Weights matrix, W=row standardized, B=binary, C=globally standardized
sids_kn1_w<-nb2listw(sids_kn1, style="W")
sids_kn1_w

sids_kn2_w<-nb2listw(sids_kn2, style="W")
sids_kn2_w

sids_kn2_b<-nb2listw(sids_kn2, style="B")
sids_kn2_b


# inverse distance weighted weights
# distance neighbours (still in degrees)
dists<-nbdists(sids_kn2, coordinates(nc_sp))

# inverse distance weighted weights
idw <- lapply(dists, function(x) 1/(x/2))
sids_idw_dist_w <- nb2listw(sids_kn2, glist = idw, style = "B")


sids74_rate=nc$SID74/nc$BIR74
# Moran's I under randomization (var(I) is normal approximation with randomization)
moranSIDS<-moran.test(nc$SID74/nc$BIR74,sids_kn2_w)
moranSIDS

# Permutations, can save each randomization I and plot distribution
moranSIDS.mc<- moran.mc(sids74_rate,sids_kn2_w, nsim=999)  
dist999<- hist(moranSIDS.mc$res,freq=TRUE,col="light blue",main="Permutation Test for Moran's I - 999 permutations",breaks=75)
lines(moranSIDS.mc$statistic,max(dist999$counts),type="h",col="red",lwd=2)

# Plot to see what distance matrix produces significant Moran's I
# all.dists<-dist(coordinates(nc))
# Increase number of groups to capture smaller distances (e.g. nbclass=50)
# Distance method not optimal for irregular areal units, better for grids

corrdist.mi<-correlog(coordinates(nc_sp),nc$SID74/nc$BIR74, method="Moran", nbclass=15)
plot(corrdist.mi,main="Moran's I for SIDS rate Correlogram, Distance Lags")

# can also use Geary's c method for correlog
corrdist.gc<-correlog(coordinates(nc_sp),nc$SID74/nc$BIR74,method="Geary", nbclass=15)
plot(corrdist.gc,main="Geary's C for SIDS rate Correlogram, Distance Lags")

library(pgirmess)
corrneigh<-sp.correlogram(sids_kn2, nc$SID74/nc$BIR74, order=15, method="I", zero.policy=TRUE)
plot(corrneigh,main="Moran's I for SIDS rate Correlogram, Neighbour Lags")


# Moving to Local Indexes of Spatial Autocorrelation
# Moran Scatterplot tells us potential outliers, divides into quadrants
mscat<-moran.plot(sids74_rate,sids_kn2_w, zero.policy=T, pch=16, col="black",cex=.5, quiet=F, labels=as.character(nc$names), main="Moran Scatterplot")

# Plot counties that are in different quadrants of High-Low SIDS and Low-High SIDS
infl<- apply(mscat$is.inf, 1, any)

# Make break points for sids rates and lagged sids rates

# Break rates into low and high
sids74_LH_breaks <- cut(sids74_rate, breaks=c(min(sids74_rate), mean(sids74_rate), max(sids74_rate)), labels=c("L", "H"), include.lowest=TRUE)
sids_lag <- lag.listw(sids_kn2_w, sids74_rate)
# Break lagged rates into low and high 
sids74_lag_breaks <- cut(sids_lag, breaks=c(min(sids_lag), mean(sids_lag), max(sids_lag)), labels=c("L", "H"), include.lowest=TRUE)
# Make quadrants LL,LH,HL,HH
lhlh <- interaction(sids74_LH_breaks, sids74_lag_breaks, infl, drop=TRUE)
# Identify low-low, low-high, etc colours
cols <- rep(1, length(lhlh))
cols[lhlh == "H.L.TRUE"] <- 2
cols[lhlh == "L.H.TRUE"] <- 3
cols[lhlh == "H.H.TRUE"] <- 4

col.palette <- colorRampPalette(c("white","orange", "red"), space = "rgb")
plot(nc_sp, col=col.palette(4)[cols], main="Morans I")
legend("bottom",horiz=TRUE , legend=c("None", "HL", "LH", "HH"), fill= col.palette(4), bty="n", cex=0.8, y.intersp=0.8)

# Local Moran's I (similar for local Geary's c)

moranLocSIDS<-localmoran(sids74_rate,sids_kn2_w,p.adjust.method="bonferroni")

moranLocSIDS2<-localmoran(nc$SID74,sids_kn2_w,p.adjust.method="bonferroni")

moranLocSIDS3<-localmoran(sids74_rate,sids_idw_dist_w,p.adjust.method="bonferroni")

# plotting the results
# Looking at the z-scores

nc_sp$locM<-(moranLocSIDS[,1])
morans_breaks <- quantile(nc_sp$locM,seq(0,1,1/5))                   
col_ramp <-col.palette(length(morans_breaks)) 

# Statistically significant local Moran I as points
LocMoranTest = moranLocSIDS[,5]<=0.05

plot(nc_sp,col=col_ramp[findInterval(moranLocSIDS[,1],morans_breaks,all.inside=TRUE)])
points(coordinates(nc_sp)[,1][LocMoranTest],coordinates(nc_sp)[,2][LocMoranTest], col="blue", cex=1, pch=18)  



#Getis-Ord G-test
# G* local Use binary weights
# recall that the Getis-Ord G statistic differs from Moran's I or Geary C 
# by the fact that it makes a difference between high-high correlation and
# low-low correlation (both of which Moran's I  treats as positive 
# autocorrelation)
# only binary non row-standardized weight matrices work for the G-test

nc$GstarSIDS<-localG(sids74_rate,sids_kn2_b,return_internals=TRUE)
plot(nc[18])

# Graphing the G-stat
gstar_breaks <- quantile(GstarSIDS,seq(0,1,1/5))                   
col_ramp <-col.palette(length(gstar_breaks))   

# picks colors corresponding for each interval out of the "cm" color palette
plot(nc_sp,border="lightgray",col=col.ramp[findInterval(GstarSIDS,gstar_breaks,all.inside=TRUE)])
# maps the shape file and fills each polygon with a color corresponding to 
# the level of G at that location
text(coordinates(nc_sp)[,1],coordinates(nc_sp)[,2], round(GstarSIDS, digits=1), cex=0.5)  
# writes the G-statistic in each polygon
legend("bottom",horiz=TRUE,legend=(round(gstar_breaks,digits=4)),fill=col_ramp,bty='n',cex=0.5)
title("Getis-Ord G statistic for SIDS rates")

