# Projecting to albers equal area optimized for California
library(maps)
library(geoR)
library(proj4)
library(ggplot2)


met_data<-read.csv("/Users/mf/Dropbox (University of Southern California)/Courses/PM569/data/harvey_825.csv")
#define the projection
projutm14<-"+proj=utm +zone=14 +north +datum=WGS84 +units=km" # TX UTM proj

#project point locations (lat,lon) to x,y 
newcoords<-project(as.matrix(cbind(met_data$lon, met_data$lat)), proj=projutm14)
met_data$x<-newcoords[,1]
met_data$y<-newcoords[,2]

#test that distances are in km
summary(dist(newcoords))

#project the state outline for mapping
#get out the map coordinates

TX <- data.frame(map("state","texas", plot=FALSE)[c("x","y")])

#project the map coordinates
newcoordsTX<-project(TX, proj=projutm14)

#map
ggplot() + geom_path(data = data.frame(newcoordsTX), aes(x=x, y = y))+
  geom_point(data = met_data, mapping = aes(x = x, y = y, color=temp))+
  coord_fixed(1.3)


# Note: can also use a geodata object 
plot(newcoordsTX,type='l')
geodata<-as.geodata(met_data,coords.col=c(13,14),data.col=6)
points(geodata,cex.max=0.8, col= rev(heat.colors(20)), pt.divide="equal")
lines(newcoordsTX)


