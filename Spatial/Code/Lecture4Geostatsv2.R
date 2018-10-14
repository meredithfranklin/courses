library(geoR) # kriging
library(fields) # plotting images
library(proj4)
library(mgcv) # smoothing


# Import data

setwd('/Users/mf/Dropbox (University of Southern California)/Courses/PM569/data')

pm_dat<-read.csv("pmdata.csv")
co2_dat<-read.csv("co2_ts.csv")

# Creation of random elevation variable to use as a spatially varying covariate in universal kriging
pm_dat$elev<-rnorm(367,1000,400)

# define the projection, distance units are km
projalbers<-"+proj=aea +lat_1=29.5 +lat_2=40.5 +lat_0=23 +lon_0=-96.0 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km"

#project point locations (lat,lon) to x,y 
newcoords<-project(as.matrix(cbind(pm_dat$lon, pm_dat$lat)), proj=projalbers)

pm_dat$x=newcoords[,1]
pm_dat$y=newcoords[,2]


# create geodata object with projected coordinates
pm_dat_geo<-as.geodata(pm_dat,coords.col=c(5,6), data.col=3, covar.col=4, covar.names="elev")
plot(pm_dat_geo)


#create grid for kriging expand.grid() in fields library
#res is the length in which we expand the max and min coordinates, 
# so we go 60 units in the x direction (from min to max x) and similarly for y
res<-60
xs<-seq(min(pm_dat$x),max(pm_dat$x),len=res)
ys<-seq(min(pm_dat$y),max(pm_dat$y),len=res)
myGrid<-expand.grid(xs,ys)
names(myGrid)<-c('x','y')

#Kriging, need the as.geodata object

# recall before, get the estimates of covariance parameters from eyeballing
vario<-variog(pm_dat_geo,uvec=seq(0,1000,l=20),option="bin",estimator.type="modulus")
plot(vario)

# profile likelihood

# Semivariogram fitting
# Note: when attempted with projected coordinates, had a lot of trouble finding the range. 
# Gaussian seems to fit best (matern not a reasonable range).
mlfit_gau<-likfit(pm_dat_geo,ini.cov.pars=c(20,400), nugget=5,kappa=1,cov.model='gaussian',fix.kappa=FALSE, fix.nugget=FALSE, lik.method='REML')
summary(mlfit_gau)
lines(mlfit_gau,col='red')

# profile method
ini.par <- cbind(seq(10, 30, l = 10), seq(300,700, l=100))

# this may take a long time
mlfit_gau2<-likfit(pm_dat_geo,ini.cov.pars=ini.par, cov.model='gaussian', fix.nugget=FALSE, lik.method='REML')
summary(mlfit_gau2)


# this too may take a long time
prof <- proflik(mlfit_gau, geodata = pm_dat_geo, sill.values = seq(20,30,l=5), range.values = seq(200, 300, l=5), nugget.values=seq(1,10,l=5),
                uni.only = FALSE) 

# Kriging by Maximum likelihood
# Simple kriging (constant known mean must be specified)
KCsimple<-krige.control(type.krige="sk",obj.m=mlfit_gau,beta=20)

# locations= are the unobserved locations where we want to make predictions
simple_krig<-krige.conv(pm_dat_geo,locations=myGrid,krige=KCsimple, output=output.control(signal=TRUE))

#plot on grid
image.plot(xs,ys,matrix(simple_krig$predict,res,res,byrow=FALSE),col=tim.colors(32), main="Simple Kriging")
#to save kriging predicted vals
myGrid2<-data.frame(myGrid,simple_krig$predict, simple_krig$krige.var)

#ordinary kriging, default in krig.control

KCord<-krige.control(type.krige='ok',obj.m=mlfit_gau)
ordinary_krige<-krige.conv(pm_dat_geo,locations=myGrid,krige=KCord,output=output.control(signal=TRUE))

#plot on grid
image.plot(xs,ys,matrix(ordinary_krige$predict,res,res,byrow=FALSE),col=tim.colors(32), main="Ordinary Kriging")
# plot standard errors
image.plot(xs,ys,matrix(sqrt(ordinary_krige$krige.var),res,res,byrow=FALSE),col=tim.colors(32))


# adding measurement locations as points to the map
points(pm_dat$x,pm_dat$y,pch=19,cex=0.5)

# add add map outline to image.plot
states<-data.frame(map("state", plot=FALSE)[c("x","y")])
#project the map coordinates
newcoordsStates<-project(states, proj=projalbers)
map(newcoordsStates,add=TRUE)



#universal kriging with trend, need to redo likfit with trend
# check trends first
summary(lm(pm~x+y,data=pm_dat))
summary(lm(pm~x+I(x^2)+y+I(y^2)+ I(x*y),data=pm_dat))

universal_mlfit_gau=likfit(pm_dat_geo,ini.cov.pars=c(20,400), nugget=5, fix.nugget=FALSE, cov.model='gaussian', lik.method='ML',trend= ~pm_dat_geo$coords[,1] )
summary(universal_mlfit_gau)

KCtrend<-krige.control(obj.m=universal_mlfit_gau,trend.d=~pm_dat_geo$coords[,1],trend.l= ~myGrid$x)
universal_krige<-krige.conv(pm_dat_geo,locations=myGrid,krige=KCtrend,output=output.control(signal=TRUE))

#plot on grid
image.plot(xs,ys,matrix(universal_krige$predict,res,res,byrow=FALSE),col=tim.colors(32))
image.plot(xs,ys,matrix(sqrt(universal_krige$krige.var),res,res,byrow=FALSE),col=tim.colors(32))

#universal kriging with covariate
# check regression first
mod_elev<-lm(pm~elev,dat=pm_dat)
summary(mod_elev)

universal_mlfit_gau2=likfit(pm_dat_geo,ini.cov.pars=c(20,400), nugget=5, fix.nugget=FALSE, cov.model='gaussian', lik.method='REML',trend= ~elev)
summary(universal_mlfit_gau2)

KCtrend2<-krige.control(obj.m=universal_mlfit_gau2,trend.d= ~elev,trend.l=~elev)
universal_krige2<-krige.conv(pm_dat_geo,locations=myGrid,krige=KCtrend2,output=output.control(signal=TRUE))
# Note we cannot do this because we don't have trend values (elevation) at all of the prediction points.

# universal krigging with trend and covariate predict at point locations
universal_mlfit_gau3=likfit(pm_dat_geo,ini.cov.pars=c(20,400), nugget=5, fix.nugget=FALSE, cov.model='gaussian', lik.method='REML',trend= ~elev)
summary(universal_mlfit_gau3)

KCtrend3<-krige.control(obj.m=universal_mlfit_gau3,trend.d= ~elev+x+y,trend.l=~elev+x+y)
universal.krige3<-krige.conv(pm_dat_geo,locations=myGrid,krige=KCtrend3,output=output.control(signal=TRUE))




# Cross validation 
# Split data into test and training sets
set.seed(123)
s<-sample(1:nrow(pm_dat),257,replace=FALSE)
temp<-1:length(pm_dat$pm)
temp[s]<-0
testset<-temp[temp!=0]
test<-pm_dat[testset,]
train<-pm_dat[s,]
# build kriging model on train, prediction locations are test
# compare the observed and predictions


####################################################################################
# Interpolation (IDW) and spatial prediction in the mean (gam)

# look at distances with two functions dist and rdist. rdist is in fields library
# distances within PM locations 
# use projected coordinates here

pm_dists<-dist(newcoords)
# distances between PM locations and grid
pm_new_dists<-rdist(newcoords,myGrid)

#inverse distance weighting function

idw<-function(data,locs,newlocs,p){
  dists<-rdist(newlocs,locs)
  return(((dists^(-p))%*%data)/((dists^(-p))%*%rep(1,length(data)))) 
}

# testing different p's
idwPred2<-idw(pm_dat$pm,cbind(pm_dat$x,pm_dat$y),myGrid,p=1)
idwPred6<-idw(pm_dat$pm,cbind(pm_dat$x,pm_dat$y),myGrid,p=6)
idwPred12<-idw(pm_dat$pm,cbind(pm_dat$x,pm_dat$y),myGrid,p=12)

# plotting the results on our grid
image.plot(xs,ys,matrix(idwPred2,res,res),col=tim.colors(32))
points(pm_dat$x,pm_dat$y,pch=19,cex=0.3)

image.plot(xs,ys,matrix(idwPred6,res,res),col=tim.colors(32))
points(pm_dat$x,pm_dat$y,pch=19,cex=0.3)

image.plot(xs,ys,matrix(idwPred12,res,res),col=tim.colors(32))
points(pm_dat$x,pm_dat$y,pch=19,cex=0.3)


# Basis and splines
# CO2 data for time series illustration of basis functions
plot(co2_dat$Decimal.date,co2_dat$CO2,type='l')
co2_dat2<-co2_dat[co2_dat$Year==2005,]
co2_dat3<-co2_dat[co2_dat$Year==2007,]
plot(co2_dat3$CO2)
# Simple polynomial basis
lm_co2=lm(CO2~Month+I(Month^2)+I(Month^3),data=co2_dat2)

# Fit polynomial result
Month<-seq(1:12)
y <- 372.921212 + 5.532769 * Month -0.978039* Month^2 + 0.047024*Month^3

plot(co2_dat2$CO2, ylab="CO2, ppb",xlab="Month")
lines(y)

# Use polynomial result to predict for different year
poly_pred<-predict(lm_co2,newdat=co2_dat3)

plot(co2_dat3$CO2,ylim=c(370,390))
lines(poly_pred,col='red')
lines(y,col='blue')

# Using cubic regression spline bases with 4 knots
library(mgcv)
gam_co2<-gam(CO2~s(Month,bs="cr"),data=co2_dat2)
plot(gam_co2)
summary(gam_co2)


# Thin plate splines (2-D)
# Let gam() find number of knots
gam_pm<-gam(pm~s(x,y,bs="ts",k=100),data=pm_dat)
plot(gam_pm)
summary(gam_pm)


# Use fitted gam object to predict smooth in x,y at grid locations

pred_gam<-predict.gam(gam_pm,myGrid,se.fit=TRUE)

# Plot predictions
image.plot(xs,ys,matrix(pred_gam$fit,res,res),col=tim.colors(32))
points(pm_dat$x,pm_dat$y,pch=19,cex=0.3)

# Standard errors of predictions
x11() # gives a new plotting window, not required
image.plot(xs,ys,matrix(pred_gam$se.fit,res,res),col=tim.colors(32))


# Change number of knots so it is more or less smooth 
gam_pm_new<-gam(pm~s(x,y,bs="ts",k=250,fx=TRUE),data=pm_dat)
summary(gam_pm_new)
plot(gam_pm_new)
pred_gam_new<-predict.gam(gam_pm_new,myGrid,se.fit=TRUE)
image.plot(xs,ys,matrix(pred_gam_new$fit,res,res),col=tim.colors(32))

gam_pm_new2<-gam(pm~s(x,y,bs="ts",k=50),data=pm_dat)
summary(gam_pm_new2)
pred_gam_new2<-predict.gam(gam_pm_new2,myGrid,se.fit=TRUE)
image.plot(xs,ys,matrix(pred_gam_new2$fit,res,res),col=tim.colors(32))


# Can also add linear covariates
gam_pm_new3<-gam(pm~s(elev)+s(x,y,bs="ts",k=100),data=pm_dat)
summary(gam_pm_new3)
pred_gam_new3<-predict.gam(gam_pm_new3,myGrid,se.fit=TRUE)
image.plot(xs,ys,matrix(pred_gam_new3$fit,res,res),col=tim.colors(32))
# Note we also have the issue here where we can't predict where don't have elevation.


