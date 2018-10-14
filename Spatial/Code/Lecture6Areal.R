library(spdep) # for areal data analysis
library(pgirmess) # for correll function, Moran's I correllogram
library(maptools) # for reading arcgis shapefiles
library(foreign) # for reading dbf files (data table associated with shapefile)

#Shapefile contains all of the data and polygon information
nc_file <- system.file("etc/shapes/sids.shp", package = "spdep")[1]
getinfo.shape(nc_file)

llCRS<-CRS("+proj=longlat +datum=NAD83")

nc <- readShapeSpatial(nc_file, ID = "FIPSNO",proj4string=llCRS)


sids79.rate <- nc@data$SID79/nc@data$BIR79
# Define categories for plotting
sids79.rate.breaks <-quantile(sids79.rate, seq(0,1,1/5))

# Define colours for plotting 
col.palette <- colorRampPalette(c("white","orange", "red"), space = "rgb")
# Pick succesive colours from our palette and assign to break points
col.ramp <-col.palette(length(sids79.rate.breaks))   

# Plot data, findInterval assigns the breaks to the data 
plot(nc,border="lightgray",col=col.ramp[findInterval(sids79.rate, sids79.rate.breaks, all.inside=TRUE)])
legend("bottom",horiz=TRUE,legend=(round(sids79.rate.breaks,digits=4)),fill=col.ramp,bty='n',cex=0.5)
#legend("bottom",legend=(round(sids79.rate.breaks,digits=4)),fill=col.ramp,bty='n',cex=0.6)

nc.counties <- fortify(nc,region="CNTY_ID")
# merge data associated with each polygon by.x is the merge variable in the shapefile, by.y is the merge variable in the dataset
nc.counties.dat  = merge(nc.counties, nc@data, by.x = 'id', by.y = 'CNTY_ID')
nc.counties.dat$sids79.rate <- nc.counties.dat$SID79/nc.counties.dat$BIR79
sids79.rate.map<-qmap("north carolina", zoom = 6, maptype = 'hybrid') +  
  geom_polygon(aes(x = long, y = lat, group = group, fill=sids79.rate), data = nc.counties.dat, colour = 'white', size = .3) +
  scale_fill_gradientn(colours=c("navyblue","darkmagenta","darkorange"))
sids79.rate.map

# k nearest neighbours
# knearneigh() creates matrix with index for the regions belonging to knn
# knn2nb() creates neighbourhood list
IDs<-row.names(as(nc, "data.frame"))
IDs<-nc@data$CNTY_ID

sids.kn1<-knn2nb(knearneigh(coordinates(nc), k=1, RANN=FALSE), row.names=IDs)
sids.kn2<-knn2nb(knearneigh(coordinates(nc), k=2, RANN=FALSE), row.names=IDs)
sids.kn3<-knn2nb(knearneigh(coordinates(nc), k=3, RANN=FALSE), row.names=IDs)
sids.kn4<-knn2nb(knearneigh(coordinates(nc), k=4, RANN=FALSE), row.names=IDs)
sids.kn10<-knn2nb(knearneigh(coordinates(nc), k=10, RANN=FALSE), row.names=IDs)

# Weights matrix, W=row standardized, B=binary, C=globally standardized

sids.kn2.w<-nb2listw(sids.kn2, style="W")
sids.kn2.w
sids.kn2.b<-nb2listw(sids.kn2, style="B")
sids.kn2.b

sids.kn3.w<-nb2listw(sids.kn3, style="W")
sids.kn3.w
sids.kn3.b<-nb2listw(sids.kn3, style="B")
sids.kn3.b

# distance neighbours (still in degrees)
ndist<-unlist(nbdists(sids.kn2, coordinates(nc)))
summary(ndist)
max.dist.kn2<-max(ndist)

dists<-nbdists(sids.kn2, coordinates(nc))

# inverse distance weighted weights
idw <- lapply(dists, function(x) 1/(x/2))
sids.dist.idw.w <- nb2listw(sids.kn2, glist = idw, style = "B")

# Weights matrix, W=row standardized, B=binary, C=globally standardized
sids.dist2<-dnearneigh(coordinates(nc), d1=0, d2=1*max.dist.kn2, row.names=IDs)
sids.dist2.w<-nb2listw(sids.dist2, style="W")
sids.dist2.w


# Moran's I under randomization (var(I) is normal approximation with randomization)
moranSIDS<-moran.test(sids79.rate,sids.kn2.w)
moranSIDS

# Permutations, can save each randomization I and plot distribution
moranSIDS.mc<- moran.mc(sids79.rate,sids.kn2.w, nsim=999)  
dist999<- hist(moranSIDS.mc$res,freq=TRUE,col="light blue",main="Permutation Test for Moran's I - 999 permutations",breaks=75)
lines(moranSIDS.mc$statistic,max(dist999$counts),type="h",col="red",lwd=2)

# Plot to see what distance matrix produces significant Moran's I
# all.dists<-dist(coordinates(nc))
# Increase number of groups to capture smaller distances (e.g. nbclass=50)
# Distance method not optimal for irregular areal units, better for grids

corrdist.mi<-correlog(coordinates(nc),nc$SID79/nc$BIR79,method="Moran", nbclass=30)
plot(corrdist.mi,main="Moran's I for SIDS rate Correlogram, Distance Lags")

# Assuming queen and rook neighbours (neighbours by borders)
sids.nb.queen<-poly2nb(nc)
sids.nb.rook<-poly2nb(nc,queen=FALSE)
# Note queen and rook neighbours are usually very similar unless you have a regular grid
ndist<-unlist(nbdists(sids.nb.rook, coordinates(nc)))
summary(ndist)
ndist<-unlist(nbdists(sids.nb.queen, coordinates(nc)))
summary(ndist)

# Visualize neighbours
plot(nc)
plot(sids.nb.queen,coordinates(nc), add=TRUE,col="green")
plot(sids.nb.rook,coordinates(nc), add=TRUE,col="blue")

# Plot of Morans I by spatial lag (orders of neighbours)
# For example, second order is neighbours of neighbours
corrneigh<-sp.correlogram(sids.kn3, nc$SID79/nc$BIR79, order=8, method="I", zero.policy=TRUE)
plot(corrneigh,main="Moran's I for SIDS rate Correlogram, Neighbour Lags")

# Geary's c

gearySIDS<-geary.test(sids79.rate,sids.kn2.w)
#gearySIDS2<-geary.test(nc$SID79,sids.kn2.w)
# check for randomization version of geary's c


corrdist.gc<-correlog(coordinates(nc),nc$SID79/nc$BIR79,method="Geary", nbclass=15)
plot(corrdist.gc,main="Geary's C for SIDS rate Correlogram, Distance Lags")


# Moving to Local Indexes of Spatial Autocorrelation
# Moran Scatterplot tells us potential outliers, divides into quadrants
mscat<-moran.plot(sids79.rate,sids.kn2.w, zero.policy=T, pch=16, col="black",cex=.5, quiet=F, labels=as.character(nc$names), main="Moran Scatterplot")

# Plot counties that are in different quadrants
infl<- apply(mscat$is.inf, 1, any)

# Make break points for sids rates and lagged sids rates

# Break rates into low and high
sids79.LH.breaks <- cut(sids79.rate, breaks=c(min(sids79.rate), mean(sids79.rate), max(sids79.rate)), labels=c("L", "H"), include.lowest=TRUE)
sids.lag <- lag(sids.kn2.w, sids79.rate)
# Break lagged rates into low and high 
sids.lag.breaks <- cut(sids.lag, breaks=c(min(sids.lag), mean(sids.lag), max(sids.lag)), labels=c("L", "H"), include.lowest=TRUE)
# Make quadrants LL,LH,HL,HH
lhlh <- interaction(sids79.LH.breaks, sids.lag.breaks, infl, drop=TRUE)
# Identify low-low, low-high, etc colours
cols <- rep(1, length(lhlh))
cols[lhlh == "H.L.TRUE"] <- 2
cols[lhlh == "L.H.TRUE"] <- 3
cols[lhlh == "H.H.TRUE"] <- 4


plot(nc, col=col.palette(4)[cols], main="Morans ")
legend("bottom",horiz=TRUE , legend=c("None", "HL", "LH", "HH"), fill= col.palette(4), bty="n", cex=0.8, y.intersp=0.8)
text(coordinates(nc), label = nc$names, cex = 0.3)

 
# Local Moran's I (similar for local Geary's c)

moranLocSIDS<-localmoran(sids79.rate,sids.kn2.w,p.adjust.method="bonferroni")

moranLocSIDS2<-localmoran(nc$SID79,sids.kn2.w,p.adjust.method="bonferroni")

moranLocSIDS3<-localmoran(sids79.rate,sids.dist2.w,p.adjust.method="bonferroni")

# plotting the results
# Looking at the z-scores
LocMoranTest1 = (moranLocSIDS[,1])
nc$locMz<-LocMoranTest1


spplot(nc, zcol="locMz", col.regions=col.palette(20), main="Local Moran's I (|z| scores)", pretty=T)
LocMoranTest = moranLocSIDS[,5]<=0.05
points(coordinates(nc)[,1][LocMoranTest],coordinates(nc)[,2][LocMoranTest],col="blue",pch=19)
# check adding points to spplot


# Statistically significant Moran I as dots
plot(nc, col=col.ramp[findInterval(sids79.rate, sids79.rate.breaks, all.inside=TRUE)], main="SIDS rates with significant Local Moran's I")
legend("bottom",horiz=TRUE,legend=(round(sids79.rate.breaks,digits=4)),fill=col.ramp,bty='n',cex=0.5)
title("Getis-Ord G statistic for SIDS rates")

LocMoranTest = moranLocSIDS[,5]<=0.05
points(coordinates(nc)[,1][LocMoranTest],coordinates(nc)[,2][LocMoranTest],col="blue",pch=19)


#Getis-Ord G-test
# G* local Use binary weights
# recall that the Getis-Ord G statistic differs from Moran's I or Geary C 
# by the fact that it makes a difference between high-high correlation and
# low-low correlation (both of which Moran's I  treats as positive 
# autocorrelation)
# only binary non row-standardized weight matrices work for the G-test

GstarSIDS<-localG(sids79.rate,sids.kn2.b,zero.policy=NULL,spChk=NULL,return_internals=TRUE)
summary(GstarSIDS)

# Graphing the G-stat

gstar.breaks <- quantile(GstarSIDS,seq(0,1,1/5))                   
col.ramp <-col.palette(length(gstar.breaks))   
  
# picks colors corresponding for each interval out of the "cm" color palette
plot(nc,border="lightgray",col=col.ramp[findInterval(GstarSIDS,gstar.breaks,all.inside=TRUE)])
# maps the shape file and fills each polygon with a color corresponding to 
# the level of G at that location
text(coordinates(nc)[,1],coordinates(nc)[,2], round(GstarSIDS, digits=1), cex=0.5)  
# writes the G-statistic in each polygon
legend("bottom",horiz=TRUE,legend=(round(gstar.breaks,digits=4)),fill=col.ramp,bty='n',cex=0.5)
title("Getis-Ord G statistic for SIDS rates")

