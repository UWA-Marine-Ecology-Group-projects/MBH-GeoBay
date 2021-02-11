### Read data ####

library( rgdal)
library( sp)
library( raster)
library(rgeos)
library(sf)



# clear environment ---
rm(list=ls()) #clear memory

studyname <- "BOSS_GB"

# Set work directory----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # sets working directory to where this script is saved (DON't MOVE)
setwd(w.dir)

# Set sub directories----
m.dir <- "/home/anitag3/MBH-GeoBay/MBdesignGB-n3"
p.dir <- paste(w.dir,"plots",sep="/")
o.dir <- paste(w.dir,"outputs",sep="/")
d.dir <- paste(w.dir,"data",sep="/")
s.dir <- paste(m.dir,"SpatialData",sep="/")


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

#gb <- readOGR(paste(d.dir, "all.zones.GB.Boss.shp", sep='/'))
gb <- st_read(paste(d.dir, "all.zones.GB.Boss.shp", sep='/'))
HPZ <- gb[1,]
#SPZ <- gb[2,]
NPZ <- gb[2,]
#MUZ <- gb[4,]
gbs1 <- gb[3,]
gbs2 <- gb[4,]
plot(gbs2)

gbs <- st_union(gbs1, gbs2, by_feature = FALSE)
plot(gbs)

#gbs <- readOGR(paste(d.dir, "GB_shallow.shp", sep='/'))


# if you want to add control sites
# controls <- readOGR(paste(s.dir, "GB-controlsites-utm.shp", sep='/'))
# controls$Name
# plot(controls[5,])
# 
# HPZ.c1 <- controls[1,]
# NPZ.c1 <- controls[2,]
# HPZ.c2 <- controls[3,]
# NPZ.c2 <- controls[4,]
# HPZ.c3 <- controls[5,]


zones <- list() # make and empty list
# add the MP zones needed
zones$HPZ <- HPZ
zones$HPZ@data[,1] <- TRUE
zones$NPZ <- NPZ
zones$NPZ@data[,1] <- TRUE
zones$GBS <- gbs
zones$GBS@data[,1] <- TRUE

# add control zones
# zones$HPZ.c1 <- HPZ.c1
# zones$HPZ.c1@data[,1] <- TRUE
# zones$HPZ.c2 <- HPZ.c2
# zones$HPZ.c2@data[,1] <- TRUE
# zones$HPZ.c3 <- HPZ.c3
# zones$HPZ.c3@data[,1] <- TRUE
# zones$NPZ.c1 <- NPZ.c1
# zones$NPZ.c1@data[,1] <- TRUE
# zones$NPZ.c2 <- NPZ.c2
# zones$NPZ.c2@data[,1] <- TRUE

#combine
zones$allSurvArea <- st_union( zones$HPZ, zones$NPZ)
zones$allSurvArea <- st_union( zones$allSurvArea, zones$GBS)
# zones$allSurvArea <- union( zones$allSurvArea, zones$HPZ.c2)
# zones$allSurvArea <- union( zones$allSurvArea, zones$HPZ.c3)
# zones$allSurvArea <- union( zones$allSurvArea, zones$NPZ.c1)
# zones$allSurvArea <- union( zones$allSurvArea, zones$NPZ.c2)

# combine for coarse bathy
#zones$coarsebathyzones <- union(zones$HPZ, zones$HPZ.c1)
#zones$coarsebathyzones <- union(zones$coarsebathyzones, zones$NPZ.c1)

#intial look to see area
plot(gb)
plot( zones$allSurvArea, col='light grey', border='light grey', add=T)
plot( zones$HPZ, col='orange', border='green', add=T)
plot( zones$NPZ, col='green', border='dark green', add=T)
plot( zones$GBS, col='blue', border='blue', add=T)

# plot( zones$HPZ.c1, add=TRUE, col='light grey', border='orange')
# plot( zones$HPZ.c2, add=TRUE, col='light grey', border='orange')
# plot( zones$HPZ.c3, add=TRUE, col='light grey', border='orange')
# plot( zones$NPZ.c1, add=TRUE, col='light grey', border='green')
# plot( zones$NPZ.c2, add=TRUE, col='light grey', border='green')


#bathymetry from Nick Mortimer (2 Aug)
bathy <- raster(paste(s.dir, 'bathy-for-Boss.tif', sep ='/'))
#s2 <- terrain(bathy, 'slope', neighbors = 8)
#plot(s2)
plot(bathy)
plot(zones$allSurvArea, add=T)

b2 <- mask(bathy, zones$allSurvArea)
plot(b2)

bhpz <- mask(bathy, HPZ)
bnpz <- mask(bathy, NPZ)
bgbs <- mask(bathy, gbs)


# Get needed coarse bathy
#b3 <- crop(bathy, zones$coarsebathyzones)
#plot(b3)
#plot(zones$coarsebathyzones, add=T)
# increase res
#b4 <- raster::disaggregate(b3, fact = c(23.2, 27.7))
#b4 <- resample(b4, b2)

#b5 <- raster::merge(b2, b4)
#plot(b5)
#plot(zones$coarsebathyzones, add=T)

#writeRaster(b5, paste(s.dir, "bathy-for-Boss.tif", sep='/'))

gb_rasters <- list()
gb_rasters$bathy <- bathy
gb_rasters$GBS <- bgbs
gb_rasters$HPZ <- bhpz
gb_rasters$NPZ <- bnpz
#Nin_rasters$bathy <- raster( x="~/NESP/MonitoringTheme/Ningaloo19/data/Bathy2Aug/depth_195_50m_WGS84.tif")
#gb_rasters$bathy <- raster( x="~/MBHdesignGB/SpatialData/GB_CMR_bathy_utm.tif")
#gb_rasters$bathy <- mask( gb_rasters$bathy, zones$allSurvArea)
#TPI from Nick Mortimer (2 Aug) Gaussian Filtered to smooth out some artefacts
#slope <- terrain(gb_rasters$bathy, "slope")
#plot(slope)
#aspect <- terrain(gb_rasters$bathy, "aspect")
#plot(aspect)
#rough<- terrain(gb_rasters$bathy, "roughness")
#plot(rough)
#gb_rasters$slope <- slope
##TPI from Nick Mortimer (8 Aug) -- had some 'tirckery' (don't know what) smooth out some artefacts and get rid of others
#Nin_rasters$TPI_gf <- raster( "/home/fos085/NESP/MonitoringTheme/Ningaloo19/data/bathy8Aug/tpi_combined_cut_nesp_25m_WGS84.tif")
#Nin_rasters$TPI_gf <- projectRaster(Nin_rasters$TPI_gf, Nin_rasters$bathy)
#Nin_rasters$TPI_gf <- flip(Nin_rasters$TPI_gf, "y")
#Nin_rasters$TPI_gf <- mask( Nin_rasters$TPI_gf, zones$allSurvArea)

#plot(gb_rasters$slope)


###################################################
#### converting polygons to a common raster.

# r <- gb_rasters$bathy
# #plot( extent( r), add=TRUE)
# #survArea first
# HPZ_raster <- rasterize( zones$HPZ, y=r, field=zones$HPZ@data[,1], bkg.value=-999, fun="first")
# NPZ_raster <- rasterize( zones$NPZ, y=r, field=zones$NPZ@data[,1], bkg.value=-999, fun="first")
# HPZ.c1_raster <- rasterize( zones$HPZ.c1, y=r, field=zones$HPZ.c1@data[,1], bkg.value=-999, fun="first")
# HPZ.c2_raster <- rasterize( zones$HPZ.c2, y=r, field=zones$HPZ.c2@data[,1], bkg.value=-999, fun="first")
# HPZ.c3_raster <- rasterize( zones$HPZ.c3, y=r, field=zones$HPZ.c3@data[,1], bkg.value=-999, fun="first")
# NPZ.c1_raster <- rasterize( zones$NPZ.c1, y=r, field=zones$NPZ.c1@data[,1], bkg.value=-999, fun="first")
# NPZ.c2_raster <- rasterize( zones$NPZ.c2, y=r, field=zones$NPZ.c2@data[,1], bkg.value=-999, fun="first")

###################################
#convert and combine
tmp1 <- as.data.frame( gb_rasters$HPZ, xy=TRUE)
tmp2 <- as.data.frame( gb_rasters$NPZ, xy=TRUE)
tmp3 <- as.data.frame( gb_rasters$GBS, xy=TRUE)

# tmp4 <- as.data.frame( HPZ.c2_raster, xy=TRUE)
# tmp5 <- as.data.frame( HPZ.c3_raster, xy=TRUE)
# tmp6 <- as.data.frame( NPZ.c1_raster, xy=TRUE)
# tmp7 <- as.data.frame( NPZ.c2_raster, xy=TRUE)
tmp8 <- as.data.frame( gb_rasters$bathy, xy=TRUE)


GBDat <- cbind( tmp1, tmp2[,3])
GBDat <- cbind( GBDat, tmp3[,3])
#GBDat <- cbind( GBDat, tmp4[,3])
# GBDat <- cbind( GBDat, tmp5[,3])
# GBDat <- cbind( GBDat, tmp6[,3])
# GBDat <- cbind( GBDat, tmp7[,3])
GBDat <- cbind( GBDat, tmp8[,3])

colnames( GBDat) <- c("Eastern", "Northing", "HPZ", "NPZ", "GBS", "bathy")

setwd("/home/anitag3/MBH-GeoBay/MBdesignGB-n3/MBH_BOSS/data")
saveRDS( GBDat, file="GBData_forBOSS.RDS")
saveRDS( gb_rasters, file="GBRasters_forBOSS.RDS")
saveRDS( zones, file="GBZones_forBOSS.RDS")




rm( MUZ_raster, HPZ_raster, r, NPZ_raster, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)

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
