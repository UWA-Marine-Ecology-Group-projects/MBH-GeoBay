### Read data ####

library( rgdal)
library( sp)
library( raster)

# staff id 00093391
# student id 21933549

#########################
#read in the country boundary -- for plotting mostly
#coastLine <- readOGR( dsn="C:/Users/00093391/Dropbox/UWA/Research Associate/MBHpackage/Ningaloo19_Data/nsaasr9nnd_02211a04es_geo___ (1)")
#proj4string( coastLine) <- CRS("+init=epsg:4283")
#coastLine <- spTransform( coastLine, "+init=epsg:4326")

############################
#read in survey areas data
# zones is a list of polygons

gb <- readOGR( dsn="~/MBHdesignGB/SpatialData/GeoBay_CMR_UTM.shp")
HPZ <- gb[1,]
SPZ <- gb[2,]
NPZ <- gb[3,]
MUZ <- gb[4,]


zones <- list()
zones$HPZ <- HPZ
zones$HPZ@data[,1] <- TRUE
#read in the state parks
zones$NPZ <- NPZ
zones$NPZ@data[,1] <- TRUE
zones$SPZ <- SPZ
zones$SPZ@data[,1] <- TRUE
zones$MUZ <- MUZ
zones$MUZ@data[,1] <- TRUE

#combine
zones$allSurvArea <- union( zones$SPZ, zones$MUZ)
zones$allSurvArea <- union( zones$allSurvArea, zones$NPZ)
zones$allSurvArea <- union( zones$allSurvArea, zones$HPZ)

#intial look to see area
plot( zones$allSurvArea, col='orange', border='orange')
plot( zones$MUZ, col='orange', border='orange')
plot( zones$SPZ, add=TRUE, col='green', border='green')
plot( zones$HPZ, add=TRUE, col='blue', border='blue')
plot( zones$NPZ, add=TRUE, col='blue', border='blue')

#bathymetry from Nick Mortimer (2 Aug)
gb_rasters <- list()
#Nin_rasters$bathy <- raster( x="~/NESP/MonitoringTheme/Ningaloo19/data/Bathy2Aug/depth_195_50m_WGS84.tif")
gb_rasters$bathy <- raster( x="~/MBHdesignGB/SpatialData/GBmultib_lidar100mUTM_CMR2.tif")
#gb_rasters$bathy <- mask( gb_rasters$bathy, zones$allSurvArea)
#TPI from Nick Mortimer (2 Aug) Gaussian Filtered to smooth out some artefacts
slope <- terrain(gb_rasters$bathy, "slope")
plot(slope)
#aspect <- terrain(gb_rasters$bathy, "aspect")
#plot(aspect)
#rough<- terrain(gb_rasters$bathy, "roughness")
#plot(rough)
gb_rasters$slope <- slope
##TPI from Nick Mortimer (8 Aug) -- had some 'tirckery' (don't know what) smooth out some artefacts and get rid of others
#Nin_rasters$TPI_gf <- raster( "/home/fos085/NESP/MonitoringTheme/Ningaloo19/data/bathy8Aug/tpi_combined_cut_nesp_25m_WGS84.tif")
#Nin_rasters$TPI_gf <- projectRaster(Nin_rasters$TPI_gf, Nin_rasters$bathy)
#Nin_rasters$TPI_gf <- flip(Nin_rasters$TPI_gf, "y")
#Nin_rasters$TPI_gf <- mask( Nin_rasters$TPI_gf, zones$allSurvArea)

plot(gb_rasters$slope)


###################################################
#### converting polygons to a common raster.

r <- gb_rasters$bathy
#plot( extent( r), add=TRUE)
#survArea first
###       ###       ### this takes a while for fine res data  ###      ###       ###
MUZ_raster <- rasterize( x=zones$MUZ, y=r, field=zones$MUZ@data[,1], bkg.value=-999, fun="first")
SPZ_raster <- rasterize( zones$SPZ, y=r, field=zones$SPZ@data[,1], bkg.value=-999, fun="first")
HPZ_raster <- rasterize( zones$HPZ, y=r, field=zones$HPZ@data[,1], bkg.value=-999, fun="first")
NPZ_raster <- rasterize( zones$NPZ, y=r, field=zones$NPZ@data[,1], bkg.value=-999, fun="first")

###################################
#convert and combine
tmp1 <- as.data.frame( MUZ_raster, xy=TRUE)
tmp2 <- as.data.frame( SPZ_raster, xy=TRUE)
tmp3 <- as.data.frame( HPZ_raster, xy=TRUE)
tmp4 <- as.data.frame( NPZ_raster, xy=TRUE)
tmp5 <- as.data.frame( gb_rasters$bathy, xy=TRUE)
tmp6 <- as.data.frame( gb_rasters$slope, xy=TRUE)

GBDat <- cbind( tmp1, tmp2[,3])
GBDat <- cbind( GBDat, tmp3[,3])
GBDat <- cbind( GBDat, tmp4[,3])
GBDat <- cbind( GBDat, tmp5[,3])
GBDat <- cbind( GBDat, tmp6[,3])

colnames( GBDat) <- c("Eastern", "Northing", "MUZ", "SPZ", "HPZ", "NPZ", "BATHY", "SLOPE")

setwd("~/MBHdesignGB/Design3/")
saveRDS( GBDat, file="GBData_forDesign3.RDS")
saveRDS( gb_rasters, file="GBRasters_forDesign3.RDS")
saveRDS( zones, file="GBZones_forDesign3.RDS")

rm( MUZ_raster, SPZ_raster, HPZ_raster, r, NPZ_raster, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)

gc()


### for GB did not do the below ###

###################
#now for the previous BRUV drops in the area

#corrected data list
BRUVS0609_correct <- read.csv( "C:/Users/00093391/Dropbox/UWA/Research Associate/MBHpackage/Ningaloo19_Data/HistoricBRUVS/2006+2009_Ningaloo_habitat_20190806.csv")

BRUVS <- list()
BRUVS[["2006"]] <- readOGR( dsn="C:/Users/00093391/Dropbox/UWA/Research Associate/MBHpackage/Ningaloo19_Data/HistoricBRUVS/Global Archive BRUVS pre-2019/2006Ningaloomaxnmetadata.shp")
BRUVS[["2009"]] <- readOGR( dsn="C:/Users/00093391/Dropbox/UWA/Research Associate/MBHpackage/Ningaloo19_Data/HistoricBRUVS/Global Archive BRUVS pre-2019/2009Ningaloomaxnmetadata.shp")

pdf( "checkingLocationsOfClean.pdf")
plot( zones[["AMP"]])
plot( BRUVS[["2006"]], add=T)
plot( BRUVS[["2009"]], add=T, col='green')
with( BRUVS0609_correct, points( Longitude,Latitude, pch=20, col='red'))
plot( zones[["AMP"]], add=TRUE)
plot( zones[["northControl"]], add=TRUE, border='blue')
plot( zones[["southControl"]], add=TRUE, border='blue')
plot( zones[["IUCN2"]], add=TRUE, border='blue')
legend( "bottomright", legend=c("orig 2006", "orig 2009", "cleaned"), pch=c(3,3,20), col=c("black","green","red"))
dev.off()

BRUVS[["2006"]] <- BRUVS[["2006"]][BRUVS[["2006"]]@data$sample %in% BRUVS0609_correct$Sample,]
BRUVS[["2009"]] <- BRUVS[["2009"]][BRUVS[["2009"]]@data$sample %in% BRUVS0609_correct$Sample,]

BRUVS[["2010"]] <- readOGR( dsn="C:/Users/00093391/Dropbox/UWA/Research Associate/MBHpackage/Ningaloo19_Data/HistoricBRUVS/Global Archive BRUVS pre-2019/2010Ningaloomaxnmetadata.shp")
BRUVS[["2013"]] <- readOGR( dsn="C:/Users/00093391/Dropbox/UWA/Research Associate/MBHpackage/Ningaloo19_Data/HistoricBRUVS/Global Archive BRUVS pre-2019/2013Ningaloomaxnmetadata.shp")
BRUVS[["2014"]] <- readOGR( dsn="C:/Users/00093391/Dropbox/UWA/Research Associate/MBHpackage/Ningaloo19_Data/HistoricBRUVS/Global Archive BRUVS pre-2019/2014Ningaloomaxnmetadata.shp")
BRUVS[["2015"]] <- readOGR( dsn="C:/Users/00093391/Dropbox/UWA/Research Associate/MBHpackage/Ningaloo19_Data/HistoricBRUVS/Global Archive BRUVS pre-2019/2015Ningaloomaxnmetadata.shp")
BRUVS[["all"]] <- do.call( "rbind", BRUVS)

#indicators for which area which historial BRUV drop is in
BRUVS$all@data$IUCN2 <- BRUVS$all@data$NorthControl <- BRUVS$all@data$SouthControl <- BRUVS$all@data$AMP <- rep( FALSE, nrow( BRUVS$all@data))
BRUVS$all@data$IUCN2[which( apply( over( BRUVS$all, zones$IUCN2), 1, function(x) !all( is.na( x))))] <- TRUE
BRUVS$all@data$NorthControl[which( apply( over( BRUVS$all, zones$northControl), 1, function(x) !all( is.na( x))))] <- TRUE
BRUVS$all@data$SouthControl[which( apply( over( BRUVS$all, zones$southControl), 1, function(x) !all( is.na( x))))] <- TRUE
BRUVS$all@data$AMP[which( apply( over( BRUVS$all, zones$AMP), 1, function(x) !all( is.na( x))))] <- TRUE

#depth of each historic drop
BRUVS$all@data$bathy <- extract( Nin_rasters$bathy, BRUVS$all)
BRUVS$all@data$bathy[ is.na( BRUVS$all@data$bathy)] <- NA

#crop to just those that we care about.
referenceBRUVS <- as.data.frame( BRUVS$all)
referenceBRUVS <- referenceBRUVS[apply( referenceBRUVS[,c("AMP","SouthControl","NorthControl","IUCN2")], 1, any),]

print( table( referenceBRUVS$year))

write.csv( referenceBRUVS, file="referenceBruvs_forDesign2.csv", row.names=FALSE)
BRUVS <- SpatialPointsDataFrame( coords=referenceBRUVS[,c("longitude","latitude")], data=referenceBRUVS, 
                                 proj4string = CRS( proj4string( Nin_rasters$TPI_gf)))
saveRDS( referenceBRUVS, file="referenceBruvs_forDesign2.RDS")

rm( BRUVS, BRUVS0609_correct)
rm( list=lso()$OTHERS)
