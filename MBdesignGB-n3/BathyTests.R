library(raster)

b250 <- raster(paste("~/MBHdesignGB/Bathy/GB-SW_250mBathy.tif"))
b250
plot(b250)
proj4string(b250) # "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"


l1 <- raster("~/MBHdesignGB/Bathy/GBlidar1.tif")
l1

l2 <- raster("~/MBHdesignGB/Bathy/GBlidar2.tif")
l2
proj4string(l2) <- proj4string(l1)

l3 <- raster("~/MBHdesignGB/Bathy/GBlidar3.tif")
l3
proj4string(l3) <- proj4string(l1)

sr <- "+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
sr2 <- "+proj=utm +zone=50 +datum=WGS84 +units=m +south +no_defs +ellps=WGS84 +towgs84=0,0,0"

b250utm <- raster::projectRaster(b250, crs=sr2)
b250utm

plot(b250utm)
writeRaster(b250utm, "~/MBHdesignGB/Bathy/GB-SW250utm.tif")
plot(l1)
plot(l2, add=T)
plot(l3, add=T)

l4 <- raster::merge(l1,l2,l3)
plot(l4)

proj4string(l1) <- sr2

#crop
c1 <- raster::crop(b250utm, l4)
plot(c1)

test2 <- raster::disaggregate(c1, fact=c(23,28))

test2

# Mask
res1 <- raster::resample(test2, l4)
res1
plot(res1)

plot(l4, add=T)

merge1 <- raster::merge(res1, l4)
plot(merge1)

writeRaster(merge1, "~/MBHdesignGB/Bathy/GB10mUTM.tif")

merge2 <- raster::projectRaster(merge1, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(merge2)
writeRaster(merge2, "~/MBHdesignGB/Bathy/GB10mLATLON.tif")

############ South West Bathy ###########

ga <- raster("~/MBHdesignGB/Bathy/ga4858_grid1_MSL.tiff")
plot(ga)
ga
ga2 <- raster::aggregate(ga, fact=3.3)
plot(ga2)
ga2

ga3 <- projectRaster(ga2, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

plot(b250utm)
plot(ga2, add=T)

e <-  drawExtent()
bathy2 <- crop(bathy2, e)
plot(bathy2)

bathy3 <- raster::disaggregate(bathy2, fact=c(26,30))
bathy3
plot(bathy3)
proj4string(bathy3) <- proj4string(ga2)

l5 <- raster::extend(l4, bathy3)
plot(l5)

bathy4 <- raster::resample(bathy3, l5)
bathy4
plot(bathy4)

bathy5 <- raster::projectRaster(bathy4, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

writeRaster(bathy4, "~/MBHdesignGB/Bathy/GB-SW_10mUTM.tif")
writeRaster(bathy5, "~/MBHdesignGB/Bathy/GB-SW_10mLATLON.tif")

bathy6 <- projectRaster(bathy4, crs="+init=epsg:32750")
plot(bathy6)
writeRaster(ga2, "~/MBHdesignGB/Bathy/SW10mUTM.tif")
writeRaster(ga3, "~/MBHdesignGB/Bathy/SW10mLATLON.tif")

b <- raster("~/MBHdesignGB/Bathy/GBSWiso.tif")

b3 <- projectRaster(b, crs="+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
b3
plot(b3)
plot(b)

# create classification matrix
reclass_df <- c(-55, -50, 50,
                -50, -45, 45,
                -45, -40, 40,
                -40, -35, 35,
                -35, -30, 30,
                -30, -25, 25,
                -25, -20, 20,
                -20, -15, 15,
                -15, -10, 10,
                -10, -5, 5,
                -5, 0, 0)
reclass_df

# reshape the object into a matrix with columns and rows
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = TRUE)
# reclassify the raster using the reclass object - reclass_m
b_class <- reclassify(b3, reclass_m)

plot(b_class)

plot(b_class,
     breaks = c(55,50,45,40,35,30,25,20,15,10,5,0),
     col = terrain.colors(12))
writeRaster(b3, "GBSWutm.tif")
b3
plot(b3)
