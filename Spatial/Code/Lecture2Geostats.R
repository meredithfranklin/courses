#################################################
# Example R code for PM569 Lecture 2        	  #
# Geostatistics	1								                #
#################################################

library(ggplot2)
library(ggmap)
library(maps)
library(mapproj)
library(proj4)

library(geoR)

setwd('/Users/mf/Dropbox (University of Southern California)/Courses/PM569/data')

pm_dat<-read.csv("pmdata.csv")
summary(pm_dat)
hist(pm_dat$pm,breaks=25)

#create log pm
pm_dat$logpm<-log(pm_dat$pm)

# Using the geoR package for spatial analysis
# Create geodata objects for original and logged PM2.5 concentrations for exploratory analysis
pm_dat_geo<-as.geodata(pm_dat,coords.col=c(1,2), data.col=3)

pm_dat_geo_log<-as.geodata(pm_dat,coords.col=c(1,2),data.col=4)

# Plotting the geodata object is useful for exploratory analysis
plot(pm_dat_geo)

# saving the plot with log PM2.5
pdf("PMexploratory.pdf")
plot(pm_dat_geo)
dev.off()

# Plotting geoR objects as point data with different colour ramps 
points(pm_dat_geo,cex.max=0.8,col=gray(seq(1, 0.1, l=100)), pt.div="equal")

# adding labels, colours, gradients
points(pm_dat_geo,cex.max=0.9, ylim= c(40,45), col= colorRampPalette(c("darkblue", "white", "darkred"))(20), 
        pt.divide="equal",main="PM2.5 concentrations (ug/m3)", xlab="Longitude",ylab="Latitude")

points(pm_dat_geo,cex.max=0.9, pt.divide="quartiles", main="PM2.5 concentrations", xlab="Longitude",ylab="Latitude")

points(pm_dat_geo,col= rev(heat.colors(10)), x.leg= -77,y.leg=36, pt.divide="quintiles", main="PM2.5 concentrations (ug/m3)", xlab="Longitude",ylab="Latitude")

# Add state boundaries to geoR plots
map("state",add=TRUE)

# subset geodata to look at a particular region

pm_dat_geo_ss<-subset(pm_dat_geo,lon> -75  & lat<42)
points(pm_dat_geo_ss,col= rev(heat.colors(20)), pt.div="equal", main="PM2.5 concentrations (ug/m3)", xlab="Longitude",ylab="Latitude")
map("state",add=TRUE)


# Nicer maps using ggplot2 and gmaps

states <- map_data("state")

ggplot() + geom_polygon(data = states, aes(x=long, y = lat, group = group)) +
  geom_point(data = pm_dat, mapping = aes(x = lon, y = lat, color=pm)) +
  coord_fixed(1.3)

east_coast <- subset(states, region %in% c("new york", "new jersey", "connecticut","rhode island","massachusetts"))
ggplot() +geom_polygon(data = east_coast, aes(x=long, y = lat, group = group))

# using ggmap for basemap
get_map(location='boston', zoom = 15, maptype='roadmap', source='google') %>% ggmap()+
  geom_point(aes(x=lon,y=lat,color=pm), alpha=0.9,size=4, data=pm.dat)+
  scale_color_gradient(low="blue", high="red")


### PROJECTIONS 
# projecting lat and lon to utm plane
# project() is from proj4 library

projutm18<-"+proj=utm +zone=18 +north +datum=WGS84 +units=km"

newcoords<-project(as.matrix(cbind(pm_dat$lon, pm_dat$lat)), proj=projutm18)
pm_dat$x<-newcoords[,1]
pm_dat$y<-newcoords[,2]


# projecting with albers equal area
# mapproject() is from maptools library

xy<-mapproject(pm_dat$lon, pm_dat$lat, proj='albers', param=c(29.5,45.5), orientation=c(90,0,-80))
pm_dat$xx<-xy$x
pm_dat$yy<-xy$y

# Create new geodata objects with projected coordinates
pm_dat_geo_utm<-as.geodata(pm_dat,coords.col=c(5,6), data.col=3)
plot(pm_dat_geo_utm)

pm_dat_geo_albers<-as.geodata(pm_dat,coords.col=c(7,8), data.col=3)
plot(pm_dat_geo_albers)


# add state boundaries, needs to be projected
east<-map_data("state",c("Massachusetts","Maine","Kentucky","Virginia","West Virginia","Maryland","Delaware","Vermont","New Hampshire","New York",
              "Ohio","Illinois","Michigan","Wisconsin","Indiana","Pennsylvania","New Jersey"),
    projection="albers",param=c(29.5,45.5))


ggplot() + geom_polygon(data = east, aes(x=long, y = lat, group = group)) +
  geom_point(data = pm_dat, mapping = aes(x = xx, y = yy, color=pm)) +
  coord_map("albers",lat0=29.5, lat1=45.5)+
  coord_fixed(1.3)



## SEMIVARIOGRAMS
# Empirical semivariograms: cloud, bin, boxplot bin
# use UTM so distances are in meters

# cloud semivariogram (all pairwise points)
vario1<-variog(pm_dat_geo_utm,option="cloud")
vario1$max.dist

plot(vario1,xlab="Distance (h), km")

# binned semivariogram
vario2<-variog(pm_dat_geo_utm,uvec=seq(0,vario1$max.dist,l=15),option="bin")
plot(vario2,xlab="Distance (h), km")

# shorter maximum distance more bins
vario3<-variog(pm_dat_geo_utm,uvec=seq(0,25,l=20),option="bin")
plot(vario3,xlab="Distance (h), km")

# boxplots
vario4<-variog(pm_dat_geo_utm,uvec=seq(0,vario1$max.dist,l=20),bin.cloud=T)
plot(vario4,bin.cloud=T,xlab="Bin")


# Robust estimator "modulus" is preferred
vario2m<-variog(pm_dat_geo_utm,uvec=seq(0,vario1$max.dist,l=20),option="bin",estimator.type="modulus")
plot(vario2m,xlab="Distance (h), km")

vario3m<-variog(pm_dat_geo_utm,uvec=seq(0,vario1$max.dist,l=20),bin.cloud=T,estimator.type="modulus")
plot(vario3m,bin.cloud=T,xlab="Distance (h), km")


# Theoretical semivariograms by curve fitting (eyeball method)
# Exponential

plot(vario2m,xlab="Distance (h)")
curve(8+32*(1-exp(-x/1000)),from = 0 , to=vario1$max.dist, col="red",ylab="Semivariance", xlab="Distance (h)",ylim=c(0,50),add=T)

# Spherical
plot(vario2m,xlab="Distance (h)")
curve(8+32*(1.5*(x/1000)-0.5*(x/1000)^3),from = 0 , col="blue",to=vario1$max.dist, ylab="Semivariance", xlab="Distance (h)",add=T)


# Exploring anisotropy (introduction for next lecture)
vario.dir<-variog4(pm_dat_geo_utm,uvec=seq(0,vario1$max.dist,l=20),option="bin",estimator.type="modulus")
plot(vario.dir)
