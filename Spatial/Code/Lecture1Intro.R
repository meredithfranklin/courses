#################################################
# R code for PM569 Lecture 1                 	  #
# Mapping in R								                  #
#################################################
library(ggplot2)
library(ggmap)
library(tidyverse)

basemap <- get_map(location='Los Angeles', maptype="satellite", zoom = 8)
ggmap(basemap)

?get_map


# Add some locations
la_locations <- tibble(location = c("Staples Center, Los Angeles", "University of Southern California" ,"Santa Monica Pier"))
# Geocode locations
geo_la_locations <- geocode(la_locations$location)
# Make a data frame with the locations
la_locations_df <- cbind(la_locations, geo_la_locations)

# Add to map
get_map("Beverly Hills", zoom = 11, maptype="hybrid") %>% ggmap() +
  geom_point(data = la_locations_df, aes(x = lon, y = lat), color = 'red', size = 3)


# reading in external data and adding it to a map

setwd('/Users/mf/Dropbox (University of Southern California)/Courses/PM569/data')
pm_dat<-read.csv("pmdata.csv")
summary(pm_dat)
View(pm_dat)


# Convert basemap to ggplot object so we can add data
# Components of this: 1) fetch map from Google API, 2) add points from pm.dat 3) add colour bar
get_openstreetmap(location='London', zoom = 15) %>% ggmap()+
  geom_point(aes(x=lon,y=lat,color=pm), alpha=0.9,size=4, data=pm_dat)+
  scale_color_gradient(low="blue", high="red")

m<-get_map(location='Boston', zoom = 5, source="google")
ggmap(m)+
  geom_point(aes(x=lon,y=lat,color=pm), alpha=0.9,size=4, data=pm_dat)+
  scale_color_gradient(low="blue", high="red")

# play around with the zoom

# try with specific lat and lon coordinates, defining a box
coords<- c(left = -90, bottom = 35, right = -75, top = 45)
get_map(location=coords, zoom = 10, maptype='roadmap', source='google') %>% ggmap()


# go through this tutorial for more on ggmap  
# https://github.com/dkahle/ggmap
