library(raster)
library(sp)
library(rgdal)
library(plyr)
library(maptools)
library(broom)
library(ggplot2)
library(gghighlight)
library(broman)
library(extrafont)
library(ggthemes)
library(dplyr)
library(viridis)
library(RColorBrewer)


w.dir <- "~/MBHdesignGB/"
s.dir <- "~/MBHdesignGB/SpatialData/"

#### load both bathymetry files (low and hi res) ----
r <- raster(paste(s.dir, "GBmultib_lidarUTM_CMR.tif", sep='/'))
plot(r)

r2 <- raster(paste(s.dir, "GB_CMR_bathy_utm.tif", sep='/'))
plot(r2)
r2


plot(r2, col = viridis(100))
plot(r, col = magma(100), add=T)

proj4string(r)
proj4string(r2)


#### Load CMR zoning ----
gb <- readOGR(paste(s.dir, "GeoBay_CMR_UTM.shp", sep='/'))
plot(gb, add=T)
gb
gbp <- proj4string(gb)

#### load planned transects in high res ----
hiMUZ <- readOGR(paste("~/MBHdesignGB/Design3", "coords_MUZ_2020-05-18.shp", sep='/'))

hiSPZ <- readOGR(paste("~/MBHdesignGB/Design3", "coords_spz_2020-05-18.shp", sep='/'))

hiNPHP <- readOGR(paste("~/MBHdesignGB/Design3", "hires_points_2020-05-18.shp", sep='/'))
## transform to utm ----
hiMUZt <- spTransform(hiMUZ, CRS(gbp))
hiSPZt <- spTransform(hiSPZ, CRS(gbp))
hiNPHPt <- spTransform(hiNPHP, CRS(gbp))
#### Load planned transects in low res ----

lowres <- readOGR(paste("~/MBHdesignGB/Design3", "lowres_points_2020-05-18.shp", sep='/'))
lowrest <- spTransform(lowres, CRS(gbp))


## save them in a new folder ----

rgdal::writeOGR(hiMUZt, dsn="~/MBHdesignGB/Design3/other", layer="coords_hiMUZ_utm", driver="ESRI Shapefile")

rgdal::writeOGR(hiSPZt, dsn="~/MBHdesignGB/Design3/other", layer= "coords_hiSPZ_utm", driver="ESRI Shapefile", overwrite_layer=TRUE)

writeOGR(hiNPHPt , dsn="~/MBHdesignGB/Design3/other", layer=paste( "coords_hiHPNP_utm"), driver="ESRI Shapefile", overwrite_layer=TRUE)

writeOGR(lowrest , dsn="~/MBHdesignGB/Design3/other", layer=paste( "coords_lowres_utm"), driver="ESRI Shapefile", overwrite_layer=TRUE) 

writeSpatialShape(hiMUZt,"test")

#### Plot ####
plot(r)
plot(gb, add=T)
plot(hiMUZt, add=T)
plot(hiSPZt, add=T)
plot(hiNPHPt, add=T)

plot(r2, col=viridis(100))
plot(gb, add=T)
plot(lowrest, pch=20, col="black", add=T)


### Plot the sampled transects ----

muz <- readOGR(paste("~/MBHdesignGB/Design3/other/sampled_transects", "Sampled_transects_hiresMUZ.shp", sep='/'))

spz <- readOGR(paste("~/MBHdesignGB/Design3/other/sampled_transects", "Sampled_transects_hiresSPZ.shp", sep='/'))

low <- readOGR(paste("~/MBHdesignGB/Design3/other/sampled_transects", "Sampled_transects_lowres.shp", sep='/'))


#### Plot ####
pal1 <- brewer.pal(n = 20, name = "YlGnBu")
plot(r2, col=viridis(100), legend = F)
plot(r2, legend.only=TRUE, col = viridis(100), horizontal=TRUE, legend.args=list(text ="Low resolution bathymetry (m)"))
plot(r, col= terrain.colors(100), add=T)
plot(gb, add=T)
plot(muz, add=T)
plot(spz, add=T)
plot(hiNPHPt, add=T)
plot(low, add=T)
