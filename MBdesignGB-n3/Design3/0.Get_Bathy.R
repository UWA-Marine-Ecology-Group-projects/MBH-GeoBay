#### prepare bathy data ####

library( rgdal)
library( sp)
library( raster)
library( MBHdesign)

w.dir <- "~/MBHdesignGB/Design3"
s.dir <- "~/MBHdesignGB/SpatialData"

## load bathy ###
mb <- raster(paste(s.dir, "GBmultibUTM.tif", sep='/'))
plot(mb)
mb

li <- raster(paste(s.dir, "GBall_lidarUTM.tif", sep='/'))
plot(li)
li

mblarge <- raster::aggregate(mb, fact=c(4.3, 3.59))
mblarge


mbres <- raster::resample(mblarge, li)
mbres

mbl <- raster::merge(mbres, li)
plot(mbl)

writeRaster(mbl, paste(s.dir, "GBmultib_lidarUTM.tif"))

### read cmr shapefile ##

gb <- readOGR(paste(s.dir, "GeoBay_CMR_UTM.shp", sep='/'))
plot(gb)
gb

mbl2 <- projectRaster(mbl, crs="+proj=utm +zone=50 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(mbl2)
mbl2

mblcmr <- mask(mbl2, gb)
plot(mblcmr)
plot(gb, add=T)

writeRaster(mblcmr, paste(s.dir, "GBmultib_lidarUTM_CMR.tif", sep='/'))

pdf( "~/MBHdesignGB/Design3/GeoBay_MultibLidar_CMR.pdf", height=6, width=10)
par( mfrow=c(1,1), oma=c(0, 0, 3, 2))
plot(mblcmr, main = 'Geographe Bay CMR Mutiblean and Lidar bathymetry')
plot(gb, add=TRUE)
dev.off()

r <- raster(paste(s.dir, "GBmultib_lidarUTM_CMR.tif", sep='/'))
plot(r)

r2 <- raster::aggregate(r, fact=c(5,5))
plot(r2)

writeRaster(r2, paste(s.dir, "GBmultib_lidar50mUTM_CMR2.tif", sep='/'))

r3 <- raster::aggregate(r, fact=c(2,2))

writeRaster(r3, paste(s.dir, "GBmultib_lidar100mUTM_CMR2.tif", sep='/'), overwrite=T)
plot(r3)


bnpz <- raster::mask(b, zones$NPZ)
writeRaster(bnpz, paste(s.dir, "NPZmultib_lidarUTM_CMR.tif", sep='/'))

bhpz <- raster::mask(b, zones$HPZ)
writeRaster(bhpz, paste(s.dir, "HPZmultib_lidarUTM_CMR.tif", sep='/'))

bspz <- raster::mask(b, zones$SPZ)
writeRaster(bspz, paste(s.dir, "SPZmultib_lidarUTM_CMR.tif", sep='/'))

bmuz <- raster::mask(b, zones$MUZ)
writeRaster(bmuz, paste(s.dir, "MUZmultib_lidarUTM_CMR.tif", sep='/'))
