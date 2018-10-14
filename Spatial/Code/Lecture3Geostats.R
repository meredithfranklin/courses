library(geoR)
library(maps)
library(proj4)
library(fields)

# API key from Google for using google maps in library(ggmap)
# install.packages("devtools")
# library(devtools)
# devtools::install_github("dkahle/ggmap")
# get api key from google https://developers.google.com/maps/documentation/geocoding/get-api-key
# restart R studio
# load library(ggmap)
# register_google(key = "AIzaSyCME84C5t6Al2SBbYLhNPyuxg3reYwghCc")
# this is my key -- can try with this, but can also get your own key

setwd('/Users/mf/Dropbox (University of Southern California)/Courses/PM569/data')

pm_dat<-read.csv("pmdata.csv")

# project coordinates
projutm18<-"+proj=utm +zone=18 +north +datum=WGS84 +units=km"
newcoords<-project(as.matrix(cbind(pm_dat$lon, pm_dat$lat)), proj=projutm18)
pm_dat$x<-newcoords[,1]
pm_dat$y<-newcoords[,2]

# create geodata object
pm_dat_geo<-as.geodata(pm_dat,coords.col=c(4,5), data.col=3)
plot(pm_dat_geo)

# Finding good starting values for the is done visually (can check theoretical semivariograms)
# Use variog object and variofit function

# empirical semivariogram find distances
#cloud semivariogrm (ALL pairwise distances, a lot of points!)
vario1<-variog(pm_dat_geo,option="cloud")
vario1$max.dist
plot(vario1) 

#binned semivariogram
vario2<-variog(pm_dat_geo,uvec=seq(0,1000,l=20),option="bin", estimator.type="modulus")
pdf("binned_semivariogram_pm.pdf")
plot(vario2,xlab="Distance (h), km")
dev.off()

# Eyeballed theoretical semivariogram (add curve to semivariogram plot)

curve(5+30*(1-exp(-x/600)),from = 0 , to=1000, ylab="Semivariance", xlab="Distance (h)",add=T,col="red",lwd=1.5)

# fitting semivariogram by ols
# ini.cov.pars (sigma2,phi)
vfit_ols=variofit(vario2,ini.cov.pars=c(30,600),nugget=5,fix.nugget=FALSE,cov.model='exponential', weights="equal")
#plot(vario2,xlab="Distance (h), km")
# add ols fitted line to the binned semivariogram and eyeballed curve

lines(vfit_ols,col="blue",lwd=1.5)

# summarizing the fit to get parameter estimates
summary(vfit_ols)

# fitting semivariogram by wls
vfit_wls=variofit(vario2,ini.cov.pars=c(30,600),nugget=5,fix.nugget=FALSE,cov.model='exponential',weights='cressie')
#plot(vario2,xlab="Distance (h), km")
lines(vfit_wls,col="purple",lwd=1.5)


summary(vfit_wls) # SSE 593, note parameter estimates are bad!

# change model type to gaussian
vfit_wls2=variofit(vario2,ini.cov.pars=c(30,600),nugget=5,fix.nugget=FALSE,cov.model='gaussian',weights='cressie')
#plot(vario2,xlab="Distance (h), km")
lines(vfit_wls2,col="cyan",lwd=1.5)

summary(vfit_wls2) # SSE 277, much better parameter estimates

# change around values of initial parameters, max.dist, nugget to see what happens
vfit_wls3=variofit(vario2,ini.cov.pars=c(50,1000),nugget=5,max.dist=1000,cov.model='gaussian',weights='cressie')
#plot(vario2,xlab="Distance (h), km")
lines(vfit_wls3,col="purple")
summary(vfit_wls3) 

# fitting covariance functions by maximum likelihood 
# choose method='ML' or 'REML', cov.model (exponential, gaussian, matern, spherical, etc.),
# trend, be careful of range parameter

mlfit_exp=likfit(pm_dat_geo,ini.cov.pars=c(30,600),nugget=5, fix.nugget=FALSE, cov.model='exponential', lik.method='ML')
summary(mlfit_exp)

mlfit_gau=likfit(pm_dat_geo,ini.cov.pars=c(30,600),nugget=5, fix.nugget=FALSE,cov.model='gaussian', lik.method='ML')
summary(mlfit_gau)

#plot
plot(vario2,xlab="Distance (h), km")
lines(mlfit_gau)
lines(mlfit_exp,col="red")

# restricted maximum likelihood
remlfit_gau=likfit(pm_dat_geo,ini.cov.pars=c(30,600),nugget=5, fix.nugget=FALSE,cov.model='gaussian', lik.method='REML')
plot(vario2,xlab="Distance (h), km")
lines(remlfit_gau)
summary(remlfit_gau)


# including a linear trend -- does this help?
# first check linear trend in lm (linear regression). Look at regression parameter estimates, are they stat sig?
mod<-lm(pm~y+x,dat=pm_dat)
summary(mod)

# add linear trend to likelihood fit
remlfit_gau_trend=likfit(pm_dat_geo,ini.cov.pars=c(30,600),nugget=5, fix.nugget=FALSE,cov.model='gaussian', lik.method='REML',trend='1st')
summary(remlfit_gau_trend)

plot(vario2,xlab="Distance (h), km")
lines(remlfit_gau_trend)


# kriging with likelihood fit
# create grid for kriging/interpolation

res=60
xs=seq(min(pm_dat$x),max(pm_dat$x),len=res)
ys=seq(min(pm_dat$y),max(pm_dat$y),len=res)
myGrid=expand.grid(xs,ys)
names(myGrid)=c('y','x')
# names of grid coordinates must match coordinate names in dataset that is being kriged

# kiriging with no trend (cte mean)
kriged_grid=krige.conv(pm_dat_geo,locations=myGrid,krige=krige.control(obj.model=remlfit_gau), output=output.control(signal=TRUE))

# kiriging with linear trend
kriged_grid_trend=krige.conv(pm_dat_geo,locations=myGrid,krige=krige.control(obj.model=remlfit_gau_trend,trend.d='1st',trend.l='1st'), output=output.control(signal=TRUE))

#plot on grid
image.plot(xs,ys,matrix(kriged_grid$predict,res,res,byrow=FALSE),col=tim.colors(32))
map('state',add=TRUE)# project MAP!!!
points(pm_dat$x, pm_dat$y, pch=19,cex=0.4)

image.plot(xs,ys,matrix(kriged_grid_trend$predict,res,res,byrow=FALSE),col=tim.colors(32))
map('state',add=TRUE)
points(pm_dat$lon, pm_dat$lat)
