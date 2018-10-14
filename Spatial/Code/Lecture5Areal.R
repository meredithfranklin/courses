##############################################
# R Code for Lecture on Areal Data (unit 1)
##############################################

library(spdep) # for areal data analysis
library(maptools) # for reading arcgis shapefiles
library(mapproj)
library(gpclib)
library(rgeos)
#library(RColorBrewer) #for map colouring
library(ggmap)

#gpclibPermit()

# Shapefile contains all of the data and polygon information
# the NC shapefile is a file that is part of the spdep library. You can find it here (MacOS) /Library/Frameworks/R.framework/Versions/3.1/Resources/library/spdep/etc/shapes/sids.shp
nc_file <- system.file("shapes/sids.shp", package = "spData")[1]
getinfo.shape(nc_file)

# alternatively can read directly with readShapeSpatial
nc <- readShapeSpatial(system.file("shapes/sids.shp",package = "spData")[1])

# set up the projection
projNAD <- CRS("+proj=lcc +lat_1=34.33333333333334 +lat_2=36.16666666666666 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +ellps=WGS84 +datum=NAD83 +units=m +no_defs")

# read in the shapefile with maptools readShapeSpatial. Can read in with or without pre-defined projection
nc <- readShapeSpatial(nc_file, ID = "FIPSNO", proj4string=CRS("+proj=lcc +lat_1=34.33333333333334 +lat_2=36.16666666666666 +lat_0=33.75 +lon_0=-79 +x_0=609601.22 +y_0=0 +ellps=WGS84 +datum=NAD83 +units=m"))

# alternatively read, then project
# nc2 <- readShapeSpatial(nc_file, ID = "FIPSNO")
# Can project after shapefile has been read in with  proj4string
# proj4string(nc2)<-projNAD

# plot polygons and add text labels

plot(nc, axes=TRUE)
text(coordinates(nc), label = nc$NAME, cex = 0.3)


# Using ggmap to get prettier maps and overlay shapefiles
# use fortify to convert to a data.frame for use with ggplot2/ggmap and plot
nc_counties <- fortify(nc,region="FIPSNO")
# merge data associated with each polygon by.x is the merge variable in the shapefile, by.y is the merge variable in the dataset
nc_counties_dat  <- merge(nc_counties, nc@data, by.x = 'id', by.y = 'FIPSNO')
# qmap is a wrapper for ggmap and get_map, allowing for acquisition of map from google and plotting shapefile on the map in one step
# geom_polygon plots the shapefile, here specify group as the varialbe that is to be displyed.
# alpha is transperancy, size is polygon outline size, colour is the polygon boundary colour
# scale_fill_gradient specifies the colour scale. colours can be specified as well as the breaks, guide="legend" is 
# used to show quantile breaks properly

# we need to have a variable called lat for ggmap. Since we have two lat.x and lat.y, we want 
# the lat from fortifying -- rename lat.x to lat and add to dataset used in ggmap.
nc_counties_dat$lat<-nc_counties_dat$lat.x

# creating our own breaks for plotting
# 10 evenly spaced categories for sids counts
breaks_sids74 <- pretty(nc$SID74, n = 10) 
# quantile breaks (6 groups)
breaks_sids74b<-quantile(nc$SID74, seq(0,1,1/6)) # add guide="legend" to scale_fill_gradientn here.

sids74_map<-qmap("north carolina", zoom = 6, maptype = 'hybrid') +  
  geom_polygon(aes(x = long, y = lat, group = group, fill=SID74), data = nc_counties_dat, colour = 'white', alpha = .5, size = .3) +
   scale_fill_gradientn(colours=c("navyblue","darkmagenta","darkorange"),breaks = breaks_sids74, labels=format(breaks_sids74))
sids74_map
# note above you can take out breaks and labels to allow for the ggmap default.

#save map with ggsave
ggsave(sids74_map,file="SIDS74_map.pdf")

# If you don't want to use a google map in the background, you can use spplot. 

spplot(nc, "SID74", col.regions= topo.colors(10), at=breaks_sids74, col="grey", main="SIDS74 Counts")
# look at help to find different colour scales

# Areal analysis

# Shared border neighbours use polygon poly2nb()
# Queen (default) and rook
sids_nb_queen<-poly2nb(nc,queen=TRUE)
plot(nc)
plot(sids_nb_queen,coordinates(nc), add=TRUE,col="green")

#sharing full borders
sids_nb_rook<-poly2nb(nc, queen=FALSE)
plot(nc)
plot(sids_nb_rook,coordinates(nc), add=TRUE,col="blue")

# see number of neighbours
summary(card(sids_nb_queen))
card(sids_nb_rook)

diffs<-diffnb(sids_nb_queen, sids_nb_rook)

plot(nc)
plot(diffs,coordinates(nc),add=TRUE, col="red")



# k nearest neighbours
# knearneigh() creates matrix with index for the regions belonging to knn
# knn2nb() creates neighbourhood list
sids_kn1<-knn2nb(knearneigh(coordinates(nc), k=1, RANN=FALSE))
sids_kn2<-knn2nb(knearneigh(coordinates(nc), k=2, RANN=FALSE))
sids_kn4<-knn2nb(knearneigh(coordinates(nc), k=4, RANN=FALSE))


plot(nc)
plot(sids_kn1, coordinates(nc),add=TRUE,col="blue")


plot(nc)
plot(sids_kn2, coordinates(nc),add=TRUE,col="green")

plot(nc)
plot(sids_kn4, coordinates(nc), add=T,col="black")

# difference between knn1 and knn2
diffs_knn<-diffnb(sids_kn1, sids_kn2)

plot(nc)
plot(diffs_knn, coordinates(nc),add=TRUE,col="red")


# distance neighbours (still in degrees)
ndist<-unlist(nbdists(sids_kn1, coordinates(nc)))
summary(ndist)

# distance based neighbours
# median distance of knn=1 is about 30km
# creating neighbours by epsilon=30km
# d1=smallest distance, usually 0, d2=distance you want to go from the point
median_dist_kn1<-median(ndist)
sids_dist1<-dnearneigh(coordinates(nc), d1=0, d2=median_dist_kn1)
sids_dist1<-dnearneigh(coordinates(nc), d1=0, d2=0.43)


plot(nc)
plot(sids.dist1, coordinates(nc), add=T,col="blue")


# 0.75 of maximum distance
max_dist_kn1<-max(ndist)
sids_dist2<-dnearneigh(coordinates(nc), d1=0.25*max_dist_kn1, d2=1*max_dist_kn1)

plot(nc)
plot(sids_dist2, coordinates(nc), add=T,col="green")



sids_dist3<-dnearneigh(coordinates(nc), d1=max_dist_kn1, d2=1.5*max_dist_kn1)

plot(nc)
plot(sids_dist3, coordinates(nc), add=T,col="purple")



# Weights matrix, W=row standardized, B=binary, C=globally standardized
sids_kn2_w<-nb2listw(sids_kn2, style="W")
sids_kn2_w
summary(unlist(sids_kn2_w$weights))

sids_dist1_b<-nb2listw(sids_dist1, style="B")
sids_dist1_b
summary(unlist(sids_dist1_b$weights))
sids_dist1_b$weights

sids_dist1_w<-nb2listw(sids_dist1, style="W")
sids_dist1_w
summary(unlist(sids_dist1_w$weights))
sids_dist1_w$weights

sids_dist1_mm<-nb2listw(sids_dist1, style="minmax")
sids_dist1_mm
summary(unlist(sids_dist1_mm$weights))
sids_dist1_mm$weights

sids_dist1_c<-nb2listw(sids_dist1, style="C")
sids_dist1_c
summary(unlist(sids_dist1_c$weights))
sids_dist1_c$weights

# For next week (Lecture 6)
# Moran's I
moranSIDS<-moran.test(nc$SID74,sids_kn2_w)

moranSIDS2<-moran.test(nc$SID74,sids_dist1_w)

# Geary's c

gearySIDS<-geary.test(nc$SID79,sids_kn2_w)

gearySIDS2<-geary.test(nc$SID79,sids_dist1_w)

# Local Moran's I and Getis-Ord G*

moranLocSIDS<-localmoran(nc$SID79,sids_kn2_w)

gearyLocSIDS<-localG(nc$SID79,sids_kn2_w)

