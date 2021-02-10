### Create inclusion probabilities ####

## for multibeam and lidar #######
#install.packages("MBHdesign")
library( MBHdesign)
library( parallel)
library( class)
library( fields)
#install.packages("pdist")
library( pdist)
library( raster)
library( rgdal)
library( sp)
library( rgeos)

# clear environment ----
rm(list = ls())

# Set working directory ####
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
# Set sub directories----
m.dir <- "/home/anitag3/MBH-GeoBay/MBdesignGB-n3"
p.dir <- paste(w.dir,"plots",sep="/")
o.dir <- paste(w.dir,"outputs",sep="/")
d.dir <- paste(w.dir,"data",sep="/")
s.dir <- paste(m.dir,"SpatialData",sep="/")


#read in data -----

# read in the inclusion probs

inclProbs <- raster(paste(d.dir, "inclProbs_forBOSS_3.tif", sep='/'))
plot(inclProbs)
inclProbs <- setValues( inclProbs, values( inclProbs) / sum( values( inclProbs), na.rm=TRUE))
plot(inclProbs)
rootInclProbs <- inclProbs
rootInclProbs <- setValues( rootInclProbs, sqrt( values( rootInclProbs)))
zones <- readRDS(paste(d.dir, "GBZones_forBOSS.RDS", sep='/')) # this one in different folder
#Deans <- readRDS( "DeansPoints_forinNOutMP-d3.RDS")
swrast <- readRDS(paste(d.dir,"GBRasters_forBOSS.RDS", sep='/'))
#if( class( BRUVS) != "SpatialPointsDataFrame")
#Deans <- SpatialPointsDataFrame( coords=Deans[,c("Longitude","Latitude")], data=Deans, proj4string = CRS( proj4string( zones[[1]])))
#proj4string(Deans) <- proj4string(swrast$bathy)
straw.nums <- readRDS(paste(d.dir, "StrawmanNumbers_Zones.RDS", sep='/'))
names(straw.nums) <- c("HPZ", "NPZ", "GBS")
straw.nums[1]

# To change straw.nums ----
# straw.nums <- c(70, 70, 200)  #numbers of drops in and out
# straw.props <- straw.nums / sum( straw.nums) # 0.625 0.375
# names( straw.nums) <- names( straw.props) <- c("HPZ", "NPZ", "GBS")
#saveRDS( straw.nums, file="StrawmanNumbers_Zones.RDS")

################################
####  choose reference sites
####  This is a one-step sample and hence
####  uses orignal inclProbs
################################

# working out density dependent probs of inclusion
#tmpD <- spTransform( Deans, CRS="+init=epsg:3577")
#tmpDist <- as.matrix( dist( coordinates( tmpD)))
#tmpDist1 <- apply( tmpDist, 1, function(x) sum( x<1000))
#Deans@data$sampleProbs <- 1 / sqrt( tmpDist1)  #to alter the importance a little bit
#Deans@data$sampleProbs <- Deans@data$sampleProbs / sum( Deans@data$sampleProbs, na.rm=TRUE)
#Deans@data$inclProbs <- raster::extract(inclProbs, Deans)
#Deans@data$sampleProbs <- Deans@data$sampleProbs * Deans@data$inclProbs #TPI inclusion probs are zero in most places


# numRef <- rep( NA, 2)
# names( numRef) <- c("InsideMP","OutsideMP")
# 
# for( zz in c( "InsideMP","OutsideMP")){
#   myZone <- zones[[zz]]
#   #if( zz == "MUS")
#   #myZone = zones$AMP - zones$IUCN2
#   tmpDrop <- as.vector( as.matrix( over( Deans, myZone)))
#   length1 <- length(tmpDrop)
#   wna <- which(is.na(tmpDrop))
#   lengthna <- length(wna)
#   length2 <- length1-lengthna
#   numRef[zz] <- min( floor( straw.nums[zz]/2), length2)
#   #numRef[zz] <- min( floor( straw.nums[zz]), length2)
#   #numRef[zz] <- min( floor( straw.nums[zz]/2), sum( length2, na.rm=TRUE))
#   prob <- Deans@data[!is.na( tmpDrop),"sampleProbs"]
#   prob <- na.exclude(prob)
#   prob <- as.vector(prob)
#   length3 <- length(prob)
#   
#   Deans@data[!is.na( tmpDrop), "sampleProbs"][sample.int(length3, numRef[zz], prob, replace=TRUE)] <- TRUE
#   #Deans@data[!is.na( tmpDrop), "f"][sample.int(length3, numRef[zz], prob, replace=FALSE)] <- TRUE
#   
#   #Deans@data[!is.na( tmpDrop), "f"][sample.int(length2, numRef[zz], prob=Deans@data[!is.na( tmpDrop),"sampleProbs"], replace=FALSE)] <- TRUE
#   #Deans@data[!is.na( tmpDrop), "f"][sample.int( sum( tmpDrop, na.rm=TRUE), numRef[zz], replace=FALSE)] <- TRUE
# }
# 
# Deans@data$sampleProbs
# 
# # load legacy sites ----
# #legacySites <- readOGR(dsn="C:/Users/00093391/Dropbox/UWA/Research Associate/MBHpackage/Ningaloo19_Data/legacySites_2019-12-23.shp")
# 
# # Deans Points are trated as legacy sites --
# Deans2 <- Deans
# table( Deans@data$f) / table( Deans2@data$f)
# 
# plot( inclProbs)
# plot( Deans, pch=20, col ='red', add=TRUE)
# #points( coordinates( Deans)[Deans@data$f,], col='red')
# 
# #legacySites <- Deans@data[Deans@data$f,]
# #legacySites <- SpatialPointsDataFrame( coords=legacySites[,c("Longitude","Latitude")], data=legacySites, proj4string=CRS(proj4string(inclProbs)))
# 
# legacySites <- Deans # 35 legacy sites



############################
####  Spatial sample of new sites
####  from altered incl. probs.
############################

### Here use quasiSamp to get random points ####
## these points will be the center of buffer for transects ###

####  Set the seed for reproducability
#set.seed( 777)
#### HAVE NOT BEEN ABLE TO MAKE THIS FUNCTION WORK ----
newSites <- list(NPZ=NULL, HPZ=NULL, GBS=NULL)
for( zz in c("HPZ", "NPZ", "GBS")){
  print( zz)
  #the number of samples to take (specified minus the legacy number)
  numby <- floor( (straw.nums[zz]))
  #numby <- floor( (straw.nums[zz] - numRef[zz])/2)
  #numby <- floor( (straw.nums[zz] - numRef[zz]))
  #set up spatial domain
  myZone <- zones[[zz]]
  #if( zz == "AMP"){
  # myZone = zones$AMP - zones$IUCN2
  #set.seed( 747)
  #}
  #tmpIP <- mask( rootInclProbs_agg_100m, myZone)
  tmpIP <- mask( inclProbs, myZone)
  tmpIP <- crop( tmpIP, myZone)
  #take the sample of clusters based on root incl probs
  newSites[[zz]] <- quasiSamp( n=numby, potential.sites=coordinates( tmpIP), inclusion.probs=values(tmpIP), nSampsToConsider=5000)
  
  #plotting (maybe remove at a later date?)
  tmpIPFull <- mask( inclProbs, myZone)
  tmpIPFull <- crop( tmpIPFull, myZone)
  plot( tmpIPFull)
  #plot( legacySites, add=TRUE, pch=1, col='red')
  points( newSites[[zz]][,c("x","y")], pch=20, col='black')
}
newSites <- do.call( "rbind", newSites)
newSites <- SpatialPointsDataFrame( coords=newSites[,c("x","y")], data=newSites, proj4string=CRS(proj4string(inclProbs)))
#some of the spatial balance is not great...  Presumably because the balance of the reference sites is also not great...


plot(inclProbs)
points(newSites, pch = 20, cex = 0.3) # 41
plot(zones$allSurvArea, add=T)
#points(Deans, pch =20, col="red") # 42
plot(zones$allSurvArea[4,])
zones$allSurvArea@data


allz <- readOGR(paste(d.dir, "all.zones.GB.Boss.shp", sep='/'))
plot(allz, add=T)
points(newSites, pch = 20, cex = 0.3)

writeOGR(newSites, o.dir, "GB_BOSS_d3", driver = "ESRI Shapefile")




### Make sure the drops are ~ 50 m apart ----

inclProbs <- raster(paste(d.dir, "inclProbs_forBOSS.tif", sep='/'))
 
 
s <- readOGR(paste(o.dir, "GB_BOSS_d1.shp", sep='/'))
s # make sure is in UTM
s <- newSites

## calculate if 2 points fall within 25m of eachother ----
# https://gis.stackexchange.com/questions/102796/remove-points-within-x-distance

## p1 ----
p1_matrix <- gWithinDistance(s, dist = 50, byid = TRUE)
diag(p1_matrix) <- NA
p1_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p1_matrix[lower.tri(p1_matrix, diag=TRUE)] <- NA
p1_matrix

colSums(p1_matrix, na.rm=TRUE) == 0
v1 <- colSums(p1_matrix, na.rm=TRUE) == 0
s[v1, ] # 98 features left

# plot --
plot(inclProbs)
plot(s[v1, ], pch=20, col='black', cex = 0.3, add=T)


## add unique id to points ----
IDnumber <- paste0('Boss',1:250)
s$IDnumber <- IDnumber
s


writeOGR(newSites, o.dir, "GB_BOSS_d3", driver = "ESRI Shapefile", overwrite = TRUE)




# save points in utm and latlong ----
d <- raster(paste(s.dir, "GB_CMR_bathy.tif", sep='/'))
d

## Save points of design ---
#all2 <- all[v1, ]
#all2 # 70 clusters
# transform again
s2 <- spTransform(s, proj4string(d))

writeOGR(newSites, o.dir, "GB_BOSS_d3_utm", driver = "ESRI Shapefile", overwrite = TRUE)


# plot(inclProbs)
# plot(all2[all2$class=="MBH",], pch=20, col='black', add=T)
# plot(all2[all2$class=="Legacy",], pch=20, col='red', add=T)

#writeOGR(all2, d.dir, "Allclusters-filtered-d7", driver="ESRI Shapefile")


## BOSS UP TO HERE ####





## separate Deans and MBH again

MBHsites <- all2[all2$class=="MBH",] #41
Deansites <- all2[all2$class=="Legacy",] #29



###############################
####  Choose new points within clusters
####  Here I need to choose transects not points
##############################

getlocal <- function(ii){
  point <- MBHsites[ii,c("x","y")]
  r2 <- rasterize( point, inclProbs, field=1)
  pbuf <- buffer( r2, width=900) ## units are in metres
  buf <- mask( inclProbs, pbuf)
  buffer <- trim(buf, pad=0)
  return( buffer)
}


sampWithOver <- 12


fullSample <- list()
fullZones <- list()



## Get bruvs accoring to clusters ---
for( ii in 1:nrow(MBHsites)){
  tmp <- getlocal(ii)
  fullZones[[ii]] <- rownames( MBHsites@data)[ii]
  tmpm <- raster::as.matrix(tmp)
  tmpm <- t(tmpm)
  tmpdf <- as.data.frame (
    cbind (coordinates (tmp), as.numeric (tmpm)))
  colnames(tmpdf) <- c("x", "y", "inclProbs_design1")
  tmpdf <- tmpdf[ order(tmpdf$y, tmpdf$x),]  # order ascending first by northing and then by easting
  fullSample[[ii]] <- quasiSamp( n=sampWithOver, potential.sites=coordinates(tmp), inclusion.probs=values(tmp), nSampsToConsider=10000)
  #fullSample[[ii]] <- transectSamp( n=sampWithOver, potential.sites = tmpdf[,c(1,2)], 
  #potential.sites= tmpdf[,c("x","y")],
  #inclusion.probs= incprobdf[,3],
  #inclusion.probs= tmpdf[,3],
  #control=gb.control
  #constrainedSet=gb.constraints.bool
  
  plot( tmp)
  points( fullSample[[ii]]$points[,c("x","y")], pch=20, col='red')
  #plot( legacySites, add=TRUE, pch=4, col='blue')
}


fullSample <- do.call( "rbind", fullSample[1:41])
bruvs <- SpatialPointsDataFrame( coords=fullSample[,c("x","y")], data=fullSample, proj4string=CRS(proj4string(inclProbs)))

plot(inclProbs)
points(bruvs, pch =20)



## Get bruvs accoring to Dean's points ---

getlocal <- function(ii){
  point <- Deansites[ii,c("x","y")]
  r2 <- rasterize( point, inclProbs, field=1)
  pbuf <- buffer( r2, width=900) ## units are in metres
  buf <- mask( inclProbs, pbuf)
  buffer <- trim(buf, pad=0)
  return( buffer)
}


sampWithOver <- 12


fullSampled <- list()
fullZonesd <- list()


for( ii in 1:nrow(Deansites)){
  tmp <- getlocal(ii)
  fullZones[[ii]] <- rownames( Deansites@data)[ii]
  tmpm <- raster::as.matrix(tmp)
  tmpm <- t(tmpm)
  tmpdf <- as.data.frame (
    cbind (coordinates (tmp), as.numeric (tmpm)))
  colnames(tmpdf) <- c("x", "y", "inclProbs_design1")
  tmpdf <- tmpdf[ order(tmpdf$y, tmpdf$x),]  # order ascending first by northing and then by easting
  fullSampled[[ii]] <- quasiSamp( n=sampWithOver, potential.sites=coordinates(tmp), inclusion.probs=values(tmp), nSampsToConsider=10000)
  #fullSample[[ii]] <- transectSamp( n=sampWithOver, potential.sites = tmpdf[,c(1,2)], 
  #potential.sites= tmpdf[,c("x","y")],
  #inclusion.probs= incprobdf[,3],
  #inclusion.probs= tmpdf[,3],
  #control=gb.control
  #constrainedSet=gb.constraints.bool
  
  plot( tmp)
  points( fullSampled[[ii]]$points[,c("x","y")], pch=20, col='red')
  #plot( legacySites, add=TRUE, pch=4, col='blue')
}


fullSampled <- do.call( "rbind", fullSampled[1:29])
bruvsd <- SpatialPointsDataFrame( coords=fullSampled[,c("x","y")], data=fullSampled, proj4string=CRS(proj4string(inclProbs)))



plot(inclProbs)
points(bruvs, pch =20)
points(bruvsd, pch = 20, col="red")
points(newSites, pch=20, col="blue")
points(Deans3, pch=20, col = "green")
points(Deans, pch=20, col = "orange")


# join clusters : new sites and Deans ----
nsdf <- as.data.frame(MBHsites)
ddf <- as.data.frame(Deansites)
head(nsdf)
names(nsdf)
nsdf <- nsdf[,c(1:5)]
head(nsdf)
#nsdf$class <- "MBH"
#
head(ddf)
names(ddf) 
ddf <- ddf[,c(18,19,17,4,5)]
head(ddf)
names(ddf) <- c("x",  "y", "inclusion.probabilities", "ID","class" )
head(ddf)
#ddf$class <- "Legacy"
# join
clusters <- rbind(nsdf, ddf)
str(clusters) # 60

# make sp --
coordinates(clusters) <- ~x+y
points(clusters, col='orange', pch=20)


# join BRUVs: bruvs and bruvsd ----
bdf <- as.data.frame(bruvs)
bddf <- as.data.frame(bruvsd)
head(bdf)
bdf <- bdf[,c(1:4)]
bdf$class <- "MBH"
head(bddf)
bddf <- bddf[,c(1:4)]
bddf$class <- "Legacy"
#join
bruvsall <- rbind(bdf, bddf)

# make sp --
coordinates(bruvsall) <- ~x+y
points(bruvsall, col='green', pch=20)


plot(inclProbs)
plot(bruvsall[bruvsall$class=="MBH",], pch=20, col='black', add=T)
plot(bruvsall[bruvsall$class=="Legacy",], pch=20, col='red', add=T)



# To save ----

pdf(paste(w.dir, "InNOutMP", "InNOutMP-Clusters-d3.pdf", sep='/'), height=7, width=8)
plot(inclProbs, main = "Inclusion probabilities - Clusters - Design3")
plot(bruvsall[bruvsall$class=="MBH",], pch=20, col='black', add=T)
plot(bruvsall[bruvsall$class=="Legacy",], pch=20, col='red', add=T)
#plot(swnp, add=T)
dev.off()




### plot nicely #####
library(dichromat)
library(RColorBrewer)
pal <- colorRampPalette(c("red","blue"))
pdf(paste(w.dir, "InNOutMP", "InNOutMP-Bathy-BRUVs-Clusters-d3.pdf", sep='/'), height=7, width=8)
plot(swrast$bathy, main ="Bathymetry - Clusters - D3", col = rev(brewer.pal(20, "RdYlBu")))
#plot(tpi, main ="Clustered Stereo-BRUVs - SW", col = pal1)
plot(bruvsall[bruvsall$class=="MBH",], pch=20, col='black', add=T)
plot(bruvsall[bruvsall$class=="Legacy",], pch=20, col='green3', add=T)
plot(zones$SWMP, add=T)
#plot(swnp, add=T)
#points(wp, add = T, pch=20, cex = 0.7, col='green')
dev.off()




####  Write the shape files

d.dir 

writeOGR(bruvsall, dsn=d.dir, layer="InNOutMP-BRUVS-d6", driver="ESRI Shapefile", overwrite_layer=TRUE)
#writeOGR(fullSample3, dsn=d.dir, layer=paste( "Bruvs4.8", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)
writeOGR(all2, dsn=d.dir, layer="InNOutMP-clusters-d6", driver="ESRI Shapefile", overwrite_layer=TRUE)
