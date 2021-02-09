library( rgdal)
library( raster) 
library(broman)

#coastLine <- readOGR( dsn="C:/Users/21933549/Dropbox/UWA/Research Associate/Ningaloo19_Data/nsaasr9nnd_02211a04es_geo___ (1)")
#legacySites <- readOGR( dsn="C:/Users/21933549/Dropbox/UWA/Research Associate/Ningaloo19_Data/legacySites_2019-12-23.shp")
newSites <- readOGR( dsn="~/MBHdesignGB/outputs/newSites_2020-05-15.shp")
zones <- readRDS( "~/MBHdesignGB/outputs/GBZones_forDesign1.RDS")
GB_rasters <- readRDS( "~/MBHdesignGB/outputs/GBRasters_forDesign1.RDS")
inclProbs <- raster( "~/MBHdesignGB/outputs/inclProbs_design1.RDS")
#BRUVS <- readRDS( "C:/Users/21933549/Dropbox/UWA/Research Associate/Ningaloo19_Data/referenceBruvs_forDesign2.RDS")
#if( class( BRUVS) != "SpatialPointsDataFrame")
#  BRUVS <- SpatialPointsDataFrame( coords=BRUVS[,c("longitude","latitude")], data=BRUVS, proj4string = CRS( proj4string( zones[[1]])))


############################
#plot bathy, TPI and zones

yel1 <- brocolors("crayons")["Banana Mania"]
pin1 <- brocolors("crayons")["Melon"]
pin2 <- brocolors("crayons")["Apricot"]
blu1 <- brocolors("crayons")["Cornflower"]
gree1 <- brocolors("crayons")["Magic Mint"]

pdf( "~/MBHdesignGB/outputs/GeoBay_Area.pdf", height=6, width=10)
par( mfrow=c(1,2), oma=c(0, 0, 3, 2))
plot( zones[["SPZ"]], col=yel1, border=yel1, main="Sampling Zones")
plot( zones[["MUZ"]], col=pin2, border=pin2, add=TRUE)
plot( zones[["HPZ"]], col=blu1, border=blu1, add=TRUE)
plot( zones[["NPZ"]], col=gree1, border=gree1, add=TRUE)
plot( gb, add=TRUE)
legend( "topleft", legend=c("Special Purpose Zone", "Multiple Use Zone", "Habitat Protecton Zone", "Land"), fill=c(yel1, pin2, blu1, gree1), cex=0.6, bty='n')

#plot( zones[["MUZ"]], col=grey(0.95), border=grey(0.75), main="Bathymetry (m)")
#plot( zones[["SPZ"]], col=grey(0.9), border='cyan', add=TRUE)
#plot( zones[["HPZ"]], col=grey(0.95), border='blue', add=TRUE)
#plot( zones[["NPZ"]], col=grey(0.95), border='blue', add=TRUE)
plot( GB_rasters$bathy, main="Bathymetry (m)")
plot( gb, add=TRUE)

plot( zones[["SPZ"]], col=grey(0.95), border=grey(0.75), main="TPI")
plot( zones[["MUZ"]], col=grey(0.95), border="cyan", add=TRUE)
plot( zones[["HPZ"]], col=grey(0.95), border="blue", add=TRUE)
plot( zones[["NPZ"]], col=grey(0.95), border="blue", add=TRUE)
tmp <- GB_rasters$TPI_gf
tmp[ tmp < -0.05] <- NA
tmp[tmp > 0.05] <- NA
tmp1 <- mask( tmp, zones$allSurvArea)
plot( tmp1, add=TRUE)
plot( coastLine, add=TRUE, col='red')

mtext( side=3, "Geographe Bay Sampling Area", outer=TRUE, cex=1.75)

dev.off()



pdf( "~/MBHdesignGB/Design3/SampleDesign_TransectsNotClus2_incProbs.pdf")
#not for the Design in each zone
par( mfrow=c(2,2), mar=c(1,1,1,1), oma=c(0,0,2,2))
#zoney <- zones[["MUZ"]]-zones[["IUCN2"]]
zoney <- zones[["MUZ"]]
plot( zoney, main="Multiple Use Zone")
tmp <- mask( inclProbs, zoney)
#tmp <- mask( bathy, zoney)
tmp <- crop( tmp, zoney)
plot( tmp, add=TRUE)
#zoneID <- over( legacySites, zoney)
#plot( legacySites[as.vector( !is.na( zoneID)),], add=TRUE, col='red', pch=4)
#zoneID <- over( newSites, zoney)
muzsites <- raster::intersect(newSitesTsp, zoney)
plot( muzsites, add=TRUE, col='black', pch=20)
#plot( newSites, add=TRUE, col='black', pch=20)
#plot( newSites[as.vector( !is.na( zoneID)),], add=TRUE, col='black', pch=20)
#legend( "topleft", legend=c("New Sites"), bty='n', col=c("red"), pch=c(4,20))

zoney <- zones[["SPZ"]]
plot( zoney, main="Special Purpose Zone")
tmp <- mask( inclProbs, zoney)
#tmp <- mask( bathy, zoney)
tmp <- crop( tmp, zoney)
plot( tmp, add=TRUE)
spzsites <- raster::intersect(newSitesTsp, zoney)
plot( spzsites, add=TRUE, col='black', pch=20)
#zoneID <- over( legacySites, zoney)
#plot( legacySites[as.vector( !is.na( zoneID)),], add=TRUE, col='red', pch=4)
#zoneID <- over( newSites, zoney)
#plot( newSites[as.vector( !is.na( zoneID)),], add=TRUE, col='black', pch=20)
#legend( "topleft", legend=c("Legacy Sites","New Sites"), bty='n', col=c("red","blue"), pch=c(4,20))

zoney <- zones[["HPZ"]]
plot( zoney, main="Habitat Protection Zone")
tmp <- mask( inclProbs, zoney)
#tmp <- mask( bathy, zoney)
tmp <- crop( tmp, zoney)
plot( tmp, add=TRUE)
hpzsites <- raster::intersect(newSitesTsp, zoney)
plot( hpzsites, add=TRUE, col='black', pch=20)
#zoneID <- over( legacySites, zoney)
#plot( legacySites[as.vector( !is.na( zoneID)),], add=TRUE, col='red', pch=4)
#zoneID <- over( newSites, zoney)
#plot( newSites[as.vector( !is.na( zoneID)),], add=TRUE, col='black', pch=20)
#legend( "topleft", legend=c("Legacy Sites","New Sites"), bty='n', col=c("red","blue"), pch=c(4,20))

zoney <- zones[["NPZ"]]
plot( zoney, main="National Park Zone")
tmp <- mask( inclProbs, zoney)
#tmp <- mask( bathy, zoney)
tmp <- crop( tmp, zoney)

plot( tmp, add=TRUE)
#zoneID <- over( legacySites, zoney)
#plot( legacySites[as.vector( !is.na( zoneID)),], add=TRUE, col='red', pch=4)
#zoneID <- over( newSites, zoney)
#plot( newSites[as.vector( !is.na( zoneID)),], add=TRUE, col='black', pch=20)
#legend( "topleft", legend=c("Legacy Sites","New Sites"), bty='n', col=c("red","blue"), pch=c(4,20))
npzsites <- raster::intersect(newSitesTsp, zoney)
plot( npzsites, add=TRUE, col='black', pch=20)

mtext( "Sample Locations", outer=TRUE, side=3, cex=1.75)

dev.off()


pdf( "~/MBHdesignGB/outputs/newdesign/SampleDesign_notClustTransects2.pdf", height=7, width=8)
plot( gb_rasters$bathy, main = "Geographe Bay not clustered Transects (1km)")
#plot( zones[["MUZ"]], main="All Zones", col=yel1, border=grey(0.75), add=TRUE)
#plot( zones[["SPZ"]], col=pin2, border=grey(0.75), add=TRUE)
#plot( zones[["HPZ"]], col=blu1, border=grey(0.75), add=TRUE)
#plot( zones[["NPZ"]], col=gree1, border=grey(0.75), add=TRUE)
#plot( inclProbs, add=TRUE)
#plot( gb_rasters$bathy, add=TRUE)
#plot( legacySites, add=TRUE, col='red', pch=4)
plot( newSitesTsp, add=TRUE, col='black', pch=20)
plot( gb, add=TRUE)
legend( "topleft", legend="Transects", bty='n', col="black", pch=20)
dev.off()

rm( list=lso()$OTHERS)
