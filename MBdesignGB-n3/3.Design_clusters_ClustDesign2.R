#### Design Clusters ####

library( rgdal)
library( sp)
library( raster)
library( MBHdesign)

setwd("~/MBHdesignGB/outputs/newdesign/")

########################
#read in the inclusion probs

inclProbs <- raster( x="inclProbs_clus_design2.tif")
inclProbs <- setValues( inclProbs, values( inclProbs) / sum( values( inclProbs), na.rm=TRUE))
rootInclProbs <- inclProbs
rootInclProbs <- setValues( rootInclProbs, sqrt( values( rootInclProbs)))
zones <- readRDS( "GBZones_forDesign1.RDS") # this one in different folder
#BRUVS <- readRDS( "referenceBruvs_forDesign2.RDS")
#if( class( BRUVS) != "SpatialPointsDataFrame")
#  BRUVS <- SpatialPointsDataFrame( coords=BRUVS[,c("longitude","latitude")], data=BRUVS, proj4string = CRS( proj4string( zones[[1]])))
straw.nums <- readRDS( "StrawmanNumbers_Zones.RDS")

##aggregate raster to speed up computation.  May need to revisit?
#rootInclProbs_agg <- aggregate( rootInclProbs, fact=20, fun=sum)  #0.5km between centres
#rootInclProbs_agg_100m <- aggregate( rootInclProbs, fact=4, fun=sum)  #100m between centres

################################
####  choose reference sites
####  This is a one-step sample and hence
####  uses orignal inclProbs
################################

#flag for whether to revisit site or not
#BRUVS@data$revisit <- FALSE

#working out density dependent probs of inclusion
#tmpBRUVS <- spTransform( BRUVS, CRS="+init=epsg:3577")
#tmpDist <- as.matrix( dist( coordinates( tmpBRUVS)))
#tmpDist1 <- apply( tmpDist, 1, function(x) sum( x<1000))
#BRUVS@data$sampleProbs <- 1 / sqrt( tmpDist1)  #to alter the importance a little bit
#BRUVS@data$sampleProbs <- BRUVS@data$sampleProbs / sum( BRUVS@data$sampleProbs, na.rm=TRUE)
#BRUVS@data$sampleProbs <- BRUVS@data$sampleProbs * BRUVS@data$inclProbs #TPI inclusion probs are zero in most places
#to get similar numbers of each year...
tmp <- tapply( BRUVS@data$sampleProbs, BRUVS@data$year, sum)
tmp1 <- BRUVS@data$sampleProbs
for( yy in unique( BRUVS@data$year))
  tmp1[BRUVS@data$year==yy] <- tmp1[BRUVS@data$year==yy] / tmp[yy]
tmp1[BRUVS@data$year==2006] <- 5 * tmp1[BRUVS@data$year==2006] #so 2006 sites are 5 times more likely than they should be.
tmp1 <- tmp1 / sum( tmp1)
BRUVS@data$sampleProbs <- tmp1

numRef <- rep( NA, 4)
names( numRef) <- c("NPZ", "HPZ", "SPZ", "MUZ")

####  Set the seed for reproducability
#set.seed( 727)

for( zz in c( "NPZ", "HPZ", "SPZ", "MUZ")){
  myZone <- zones[[zz]]
  if( zz == "MUS")
    myZone = zones$AMP - zones$IUCN2
  tmpDrop <- as.vector( as.matrix( over( BRUVS, myZone)))
  numRef[zz] <- min( floor( straw.nums[zz]/4), sum( tmpDrop, na.rm=TRUE))
  
  BRUVS@data[!is.na( tmpDrop), "revisit"][sample.int( sum( tmpDrop, na.rm=TRUE), numRef[zz], prob=BRUVS@data[!is.na( tmpDrop),"sampleProbs"], replace=FALSE)] <- TRUE
  #  BRUVS@data[!is.na( tmpDrop), "revisit"][sample.int( sum( tmpDrop, na.rm=TRUE), numRef[zz], replace=FALSE)] <- TRUE
}

# load legacy sites
legacySites <- readOGR(dsn="C:/Users/00093391/Dropbox/UWA/Research Associate/MBHpackage/Ningaloo19_Data/legacySites_2019-12-23.shp")
legacySites@data$year
table( legacySites@data$year) / table( BRUVS@data$year)

plot( rootInclProbs_agg)
plot( BRUVS, add=TRUE)
points( coordinates( BRUVS)[BRUVS@data$revisit,], col='red')

legacySites <- BRUVS@data[BRUVS@data$revisit,]
legacySites <- SpatialPointsDataFrame( coords=legacySites[,c("longitude","latitude")], data=legacySites, proj4string=CRS(proj4string(inclProbs)))
# 49 legacy sites

############################
####  Spatial sample of new sites
####  from altered incl. probs.
############################

### Here use quasiSamp to get random points ####
## these points will be the center of buffer for transects ###

####  Set the seed for reproducability
set.seed( 777)
#### HAVE NOT BEEN ABLE TO MAKE THIS FUNCTION WORK ----
newSites <- list(NPZ=NULL,HPZ=NULL,SPZ=NULL,MUZ=NULL)
for( zz in c("NPZ", "HPZ", "SPZ", "MUZ")){
  print( zz)
  #the number of samples to take (specified minus the legacy number)
  numby <- floor( (straw.nums[zz]))
  #numby <- floor( (straw.nums[zz] - numRef[zz])/4)
  #set up spatial domain
  myZone <- zones[[zz]]
  #if( zz == "AMP"){
   # myZone = zones$AMP - zones$IUCN2
    #set.seed( 747)
  #}
  #tmpIP <- mask( rootInclProbs_agg_100m, myZone)
  tmpIP <- mask( rootInclProbs, myZone)
  tmpIP <- crop( tmpIP, myZone)
  #take the sample of clusters based on root incl probs
  newSites[[zz]] <- quasiSamp( n=numby, potential.sites=coordinates( tmpIP), inclusion.probs=values(tmpIP), nSampsToConsider=5000)
  
  #plotting (maybe remove at a later date?)
  tmpIPFull <- mask( rootInclProbs, myZone)
  tmpIPFull <- crop( tmpIPFull, myZone)
  plot( tmpIPFull)
  #plot( legacySites, add=TRUE, pch=1, col='red')
  points( newSites[[zz]][,c("x","y")], pch=20, col='black')
}
newSites <- do.call( "rbind", newSites)
newSites <- SpatialPointsDataFrame( coords=newSites[,c("x","y")], data=newSites, proj4string=CRS(proj4string(inclProbs)))
#some of the spatial balance is not great...  Presumably because the balance of the reference sites is also not great...

###############################
####  Choose new points within clusters
####  Here I need to choose transects not points
##############################

getlocal <- function(ii){
  point <- newSites[ii,c("x","y")]
  r2 <- rasterize( point, rootInclProbs, field=1)
  pbuf <- buffer( r2, width=2000) ## units are in metres
  buf <- mask( rootInclProbs, pbuf)
  buffer <- trim(buf, pad=0)
  return( buffer)
}

#sampWithOver <- 6+2
sampWithOver <- 2+1

fullSample <- list()
fullZones <- list()

## I think in this funtion I need to change quasiSamp for TransectSamp
for( ii in 1:nrow( newSites)){
  tmp <- getlocal(ii)
  fullZones[[ii]] <- rownames( newSites@data)[ii]
  tmpm <- raster::as.matrix(tmp)
  tmpm <- t(tmpm)
  tmpdf <- as.data.frame (
    cbind (coordinates (tmp), as.numeric (tmpm)))
  colnames(tmpdf) <- c("x", "y", "inclProbs_design1")
  tmpdf <- tmpdf[ order(tmpdf$y, tmpdf$x),]  # order ascending first by northing and then by easting
  #fullSample[[ii]] <- quasiSamp( n=sampWithOver, potential.sites=coordinates( tmp), inclusion.probs=values( tmp), nSampsToConsider=10000)
  fullSample[[ii]] <- transectSamp( n=sampWithOver, potential.sites = tmpdf[,c(1,2)], 
                                    #potential.sites= tmpdf[,c("x","y")],
                          #inclusion.probs= incprobdf[,3],
                          inclusion.probs= tmpdf[,3],
                          control=gb.control
                          #constrainedSet=gb.constraints.bool
  )
  plot( tmp)
  points( fullSample[[ii]]$points[,c("x","y")], pch=20, col='red')
  #plot( legacySites, add=TRUE, pch=4, col='blue')
}

### UP TO HERE ####


fullSamplep <- do.call( "rbind", fullSample[19:36])
fullSamplep <- rbind(fullSample[[1]]$points, fullSample[[2]]$points, fullSample[[3]]$points, fullSample[[4]]$points, fullSample[[5]]$points,
                     fullSample[[6]]$points, fullSample[[7]]$points, fullSample[[8]]$points, fullSample[[9]]$points, fullSample[[10]]$points, 
                     fullSample[[11]]$points, fullSample[[12]]$points, fullSample[[13]]$points, fullSample[[14]]$points, fullSample[[15]]$points,
                     fullSample[[16]]$points, fullSample[[17]]$points, fullSample[[18]]$points)
#fullSamplet <- do.call( "rbind", fullSample$transects)
#fullSample2 <- bind_rows(fullSample, .id = "column_label")
fullSamplep$cluster <- rep( do.call( "c", fullZones), each=75) # 15 is number of point per transect
#tpoints <- fullSample$points[,c("x","y")]
write.csv(fullSamplep, "~/MBHdesignGB/outputs/newdesign/GB_design2_ClusTransects.csv")

samplepoints <- fullSamplep
coordinates(samplepoints) <- ~x+y

for( zz in c("NPZ", "HPZ", "SPZ", "MUZ")){
  plot( zones[[zz]])
  plot( inclProbs, add=TRUE)
  plot( zones[[zz]], add=TRUE)
  points( samplepoints, pch=20, col='red')
  
  #plot( legacySites, add=TRUE, pch=4, col='blue')
}

fullSample$ID <- paste( fullSample$cluster, rep( paste0( "shot.",1:6), each=nrow( newSites)), sep="_")
fullSamplesp <- SpatialPointsDataFrame( coords=fullSamplep[,c("x","y")], data=fullSamplep, proj4string=CRS(proj4string(inclProbs)))

##################################
####  Write the shape files

writeOGR( fullSamplesp, dsn="~/MBHdesignGB/outputs/newdesign", layer=paste( "ClustTransects_design2", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)
#writeOGR( legacySites, dsn="C:/Users/21933549/Dropbox/UWA/Research Associate/Ningaloo19_Data", layer=paste( "legacySites", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)

##################################
####  Changing one version of the decimal degrees to degrees, minutes, seconds

convert2dms <- function(xx, fractionsOfSeconds=1){
  #this is a bit of a hack (with the signage being dropped)
  neg.pattern <- (xx < 0)
  xx[neg.pattern] <- -xx[neg.pattern]
  
  degrees <- floor( xx) #these are all positive now
  degree.remainder <- xx %% 1
  x <- degree.remainder * 60
  minutes <- floor( x)
  minute.remainder <- x %% 1
  x <- minute.remainder * 60
  seconds <- round( x, fractionsOfSeconds) * 10^fractionsOfSeconds
  
  tmp <- paste( degrees, minutes, seconds, sep=".")
  #  tmp[neg.pattern] <- paste0("-",tmp[neg.pattern])
  
  return( tmp)
}

newSitesLocations_dms <- apply( coordinates( fullSample), 2, convert2dms)
legacySitesLocations_dms <- apply( coordinates( legacySites), 2, convert2dms)
colnames( newSitesLocations_dms) <- colnames( legacySitesLocations_dms) <- c("longitude","latitude")

newSitesLocations_dms[,"longitude"] <- paste0( newSitesLocations_dms[,"longitude"], "E")
newSitesLocations_dms[,"latitude"] <- paste0( newSitesLocations_dms[,"latitude"], "S")

legacySitesLocations_dms[,"longitude"] <- paste0( legacySitesLocations_dms[,"longitude"], "E")
legacySitesLocations_dms[,"latitude"] <- paste0( legacySitesLocations_dms[,"latitude"], "S")

#write the locations

write.table( newSitesLocations_dms, file=paste0( paste("./Designs/newSitesLocations_dms", Sys.Date(), sep="_"), ".txt"), sep=" ", row.names=FALSE, quote=FALSE)
write.table( legacySitesLocations_dms, file=paste0( paste("./Designs/legacySitesLocations_dms", Sys.Date(), sep="_"), ".txt"), sep=" ", row.names=FALSE, quote=FALSE)

rm( list=lso()$OTHERS)
rm( getlocal, convert2dms)
