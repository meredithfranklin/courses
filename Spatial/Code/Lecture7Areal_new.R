library(classInt)
library(spdep)


setwd("/Users/mf/Dropbox (University of Southern California)/Courses/PM569/data")
load("Election.RData")

# shapefile election already saved as R shapefile

data <- election
names(data)

# use queen weight matrix
elec_nb <- poly2nb(data, queen=T)
elec_mat_w <- nb2listw(elec_nb, style="W", zero.policy=TRUE)
elec_mat_b <- nb2listw(elec_nb, style="B", zero.policy=TRUE)

#data_sf <- st_as_sf(data)

# OLS Linear Model
mod_lm <- lm(Bush_pct ~ pcincome + urbrural + UNDER18 + Obese + Noins, data=data)
summary(mod_lm)

# Plot residuals

res_lm <- mod_lm$residuals
res_palette <- colorRampPalette(c("red","orange","white", "lightgreen","green"), space = "rgb")
pal <- res_palette(5)

classes <- classIntervals(res_lm, n=5, style="fixed", fixedBreaks=c(-61,-25,-5,5,25,50), rtimes = 1)
cols <- findColours(classes,pal)


plot(data,col=cols,border="grey",main="Residuals from OLS Model")
legend(x="bottom",cex=1,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),ncol=5)


# OLS Residual Autocorrelation global Moran's

moran.test(res_lm, listw=elec_mat_w, zero.policy=T)


# SAR model
mod_sar1 <- spautolm(Bush_pct ~ pcincome, data = data, listw=elec_mat_w, zero.policy=T, tol.solve=1e-12,family="SAR")
summary(mod_sar1)
res_sar1 <- residuals(mod_sar1)
# 1 time Moran's test
moran_res1<-moran.test(res_sar1, listw=elec_mat_w, zero.policy=T)

# MCMC Moran's test giving us 1000 Moran's I statistics (better when # polygons small)
moran_mc_res1<- moran.mc(res_sar1,elec_mat_w, zero.policy=T,nsim=999) 
dist999<- hist(moran_mc_res1$res,freq=TRUE,col="light blue",main="Permutation Test for Moran's I - 999 permutations",breaks=75)
lines(moran_mc_res1$statistic,max(dist999$counts),type="h",col="red",lwd=2)

classes <- classIntervals(res_sar1, n=5, style="fixed", fixedBreaks=c(-50,-25,-5,5,25,50), rtimes = 1)
cols <- findColours(classes,pal)
plot(data,col=cols, border="grey",main="Residuals from SAR (error) Model")
legend(x="bottom",cex=1,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),ncol=5)


# SAR error model - same as SAR in spautolm above, applies sp weights to y and x
mod_sar2 <- errorsarlm(Bush_pct ~ pcincome, data = data, listw=elec_mat_w, zero.policy=T, tol.solve=1e-12)
summary(mod_sar2)


# SAR lag model different function - only applies sp weights to y
mod_sar3 = lagsarlm(Bush_pct ~ pcincome, data = data, elec_mat_w, zero.policy=T, tol.solve=1.0e-30)
summary(mod_sar3)
res_sar3<-mod_sar3$residuals
moran.test(res_sar3, listw=elec_mat_w, zero.policy=T)

classes <- classIntervals(res_sar3, n=5, style="fixed", fixedBreaks=c(-50,-25,-5,5,25,50), rtimes = 1)
cols <- findColours(classes,pal)
plot(data,col=cols, border="grey",main="Residuals from SAR lag Model")
legend(x="bottom",cex=1,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),ncol=5)



# CAR model - requires binary weight matrix
mod_car <- spautolm(Bush_pct ~ pcincome, data = data, listw=elec_mat_b, zero.policy=T, tol.solve=1e-12,family="CAR")
summary(mod_car)
res_car <- residuals(mod_car)
# Moran's test
moran_resc<-moran.test(res_car, listw=elec_mat_w, zero.policy=T)
classes <- classIntervals(res_car, n=5, style="fixed", fixedBreaks=c(-50,-25,-5,5,25,50), rtimes = 1)
cols <- findColours(classes,pal)
plot(data,col=cols, border="grey",main="Residuals from CAR Model")
legend(x="bottom",cex=1,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),ncol=5)


# Random effects model
library(nlme)

# election data -- this takes a long time to run!
data$lat<-coordinates(data)[,2]
data$lon<-coordinates(data)[,1]
spatial<-corSpatial(1,form=~lon+lat,type="gaussian") # like choosing a gaussian covariance function
scor<-Initialize(spatial, as(data, "data.frame")[,c("lon","lat")], nugget=FALSE)
mod.mixed<-lme(Bush_pct ~ pcincome, random= ~1|FIPS_num, data=as(data,"data.frame"), correlation=scor, method="ML")


# SIDS data
library(maptools)
library(geoR)


nc <- readShapeSpatial(system.file("shapes/sids.shp", package = "spData")[1])

nc$lat<-coordinates(nc)[,2]
nc$lon<-coordinates(nc)[,1]

sids_rate = (nc$SID79)/nc$BIR79
nwbirth_rate = (nc$NWBIR79)/nc$BIR79

spatial<-corSpatial(1,form=~lon+lat,type="gaussian") # like choosing a gaussian covariance function
scor<-Initialize(spatial, as(nc, "data.frame")[,c("lon","lat")], nugget=FALSE)
mod_mixed<-lme(sids_rate ~ nwbirth_rate, random= ~1|CNTY_ID, data=as(nc,"data.frame"), correlation=scor, method="ML",
               control=lmeControl(singular.ok=TRUE, returnObject = TRUE))
summary(mod_mixed)
res_mixed<-mod_mixed$residuals
nc$res_mixed<-res_mixed

sids_nb <- poly2nb(nc, queen=T)
sids_mat_w <- nb2listw(sids_nb, style="W", zero.policy=TRUE)
moran_resc<-moran.test(nc$res_mixed, listw=sids_mat_w, zero.policy=T)
# not working

# compare SIDS mixed effects with SAR model
sids_nb <- poly2nb(nc, queen=T)
sids_mat_w <- nb2listw(sids_nb, style="W", zero.policy=TRUE)
mod_sar_nc <- spautolm(sids_rate ~ nwbirth_rate, data = nc, listw=sids_mat_w, zero.policy=T, tol.solve=1e-12,family="SAR")
summary(mod_sar_nc)

# Geographically weighted regression
library(spgwr)

# GWR of election data - takes a long time to run!
#bwG <- gwr.sel(Bush_pct ~ pcincome, data = data, gweight=gwr.Gauss, verbose=F)
#mod.gwr <- gwr(Bush_pct ~ pcincome, data = data, bandwidth=bwG, gweight=gwr.Gauss)
#mod.gwr <- gwr(Bush_pct ~ pcincome, data = data, bandwidth=0.2, gweight=gwr.Gauss,longlat=TRUE)

# GWR SIDS data, use rates
sids_rate = 1000*(nc$SID79+1)/nc$BIR79
hist(sids_rate)

# Freeman-Tukey transformation of the response to make normal distribution
sids_ft = sqrt(1000)*(sqrt(nc$SID79/nc$BIR79) + sqrt((nc$SID79+1)/nc$BIR79))
hist(sids_ft)

# the distribution of the transformed values is more symmetric and there is only one outlier.
# exploration of the covariate nwbirth = non-white births
nwbirth_rate = (nc$NWBIR79)/nc$BIR79
hist(nwbirth_rate)

# Freeman-Tukey transformation of the covariate to make normal distribution
nwbirth_ft = sqrt(1000)*(sqrt(nc$NWBIR79/nc$BIR79) +sqrt((nc$NWBIR79+1)/nc$BIR79))
hist(nwbirth_ft)

# GWR takes two steps, 1) estimate bandwidth
sids_bwG<-gwr.sel(sids_ft~nwbirth_ft, data=nc, gweight=gwr.Gauss, verbose=T, adapt=T)
# 2) fit model
mod_gwr <- gwr(sids_ft ~ nwbirth_ft, data = nc, bandwidth=sids_bwG, gweight=gwr.Gauss, hatmatrix = TRUE)
# if you have a pre-defined bandwidth, you only need to do step 2.
mod_gwr
# Localized Coefficients from GWR

coef <- mod_gwr$SDF$nwbirth_ft

classes_fx <- classIntervals(coef, n=5, style="fixed", fixedBreaks=seq(-0.001,0.04,by=0.005), rtimes = 1)
cols <- findColours(classes_fx,pal)

plot(nc,col=cols, border="grey",pretty=T)
legend(x="bottom",cex=1,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Local Coefficient Estimates (non-white births)",ncol=3)

# plot local R^2
localR2<-mod_gwr$SDF$localR2

classes_fx <- classIntervals(localR2, n=5, style="fixed", fixedBreaks=seq(0.047,0.31,by=0.05), rtimes = 1)
cols <- findColours(classes_fx,pal)

plot(nc,col=cols, border="grey")
legend(x="bottom",cex=1,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Local R-sq",ncol=3)


# Residuals of GWR to see if there is residual spatial autocorrelation

res_gwr <- mod_gwr$SDF$gwr.e

classes <- classIntervals(res_gwr, n=5, style="fixed", fixedBreaks=c(-2,-0.5,0,0.5,1,2), rtimes = 1)
cols <- findColours(classes,pal)

par(mar=rep(0,4))
plot(nc,col=cols, border="grey",pretty=T)
legend(x="bottom",cex=1,fill=attr(cols,"palette"),bty="n",legend=names(attr(cols, "table")),title="Residuals from GWR Model",ncol=5)

# Test residual autocorrelation with Moran's


moran.test(res_gwr, listw=sids_mat_w)



         
