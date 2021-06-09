# Get inclusion probabilities ---

## Create inclusion probabilities ####

library( rgdal)
library( sp)
library( raster)

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



#read in data

SWDat <- readRDS(paste(d.dir, "GBData_forBOSS.RDS", sep='/'))
#sw_rasters <- raster(paste(d.dir, "bathy-for-boss.tif", sep='/'))
sw_rasters <- readRDS(paste(d.dir,"GBRasters_forBOSS.RDS", sep='/'))
zones <- readRDS(paste(d.dir, "GBZones_forBOSS.RDS", sep='/'))

# GB marine park polygon ----
#swnp <- readOGR(paste(s.dir, "SW_CMR_NP.shp", sep='/'))
hpz <- zones$HPZ
plot(hpz)

npz <- zones$NPZ
plot(npz, add=T)

gbs <- zones$GBS
plot(gbs, add=T)


####  Straw man for numbers of samples in each region ----

straw.nums <- c(50, 50, 150)  #numbers of drops in and out
straw.props <- straw.nums / sum( straw.nums) # 0.625 0.375
names( straw.nums) <- names( straw.props) <- c("HPZ", "NPZ", "GBS")
#saveRDS( straw.nums, file="StrawmanNumbers_Zones.RDS")



####  Hand-picking Bathy cut points ----
####  And their numbers of drops

#Bathy.quant <- c(0,0.8,0.9,0.925,0.95,0.975,1)
#Bathy.quant <- c(0,0.25,0.5,0.8,0.9,0.925,0.95,0.975,1)
#Bathy.quant <- c(0,0.017,0.06,0.085,0.48,0.51,0.79,0.85,1)
Bathy.quant <- c(0,0.1,0.5,0.85,1)
Bathy.cuts <- quantile( sw_rasters$bathy, Bathy.quant)#c( -Inf,0.02,0.04,0.08,0.16,Inf)
Bathy.cuts
#trying to make it so there is no hand-picking (except for the hand-picked function)
tmp <- cumsum( Bathy.quant)
#Bathy.targetNums <- rep( floor( 18/8), 4)#floor( ( tmp / sum( tmp))[-1] * 200)#rep( 40, 5)#c( 20,20,30,65,65)
#Bathy.targetNums <- rep( floor( 250/8), 4) # 
#Bathy.targetProps <-  Bathy.targetNums / sum( Bathy.targetNums) # 0.5 0.5
Bathy.targetProps <- c(0.1,0.3, 0.4, 0.2)
Bathy.targetProps.hpz <- c(0.3,0.3,0.3,0.1)
Bathy.targetProps.npz <- c(0,0.02,0.98)

# Proportion of potential sites in each zone ----

SWDat_small <- SWDat[!is.na( SWDat$bathy),]
tmp <- colSums( SWDat_small[,c("HPZ", "NPZ", "GBS")], na.rm=TRUE) 
#tmp[2] <- tmp[2] - tmp[1] # so similar amount of sites in SPZ and MUZ
tmp[1] # 7495
tmp[2] # 2676 
tmp[3]

props <- tmp / nrow( SWDat_small) # inside 0.6966145 - outside 0.2459908 
props <- props / sum( props) # inside  0.5502533 - outside 0.4497467  


###################################
####  To get cut points
###################################

catB <- cut(sw_rasters$bathy, breaks=Bathy.cuts, na.rm=TRUE)
plot(catB)
plot(zones$allSurvArea, add=T)

#plot( zones$InsideMP, add=T); plot( catB, add=TRUE); plot( zones$OutsideMP, add=TRUE)
plot( catB); plot( zones$HPZ, add=T);  plot( zones$NPZ, add=TRUE);  plot( zones$GBS, add=TRUE)
#plot(swnp, add=T)

#writeRaster(catB, paste(d.dir, 'Bathy_cuts_BOSS.tif', sep='/'), overwrite=TRUE)

plot(catB)


##################################
####  Within each zone (incl probs)
####  Weight according to straw.props
##################################


inclProbs <- catB
for( zz in c( "HPZ", "NPZ", "GBS")){
  print( zz)
  #if( zz == "NPZ")
  #zoneID <- extract( x=catB, y=zones$NPZ, cellnumbers=TRUE)
  #zoneID <- extract( x=catB, y=zones$MUZ-zones$NPZ, cellnumbers=TRUE)
  #else
  zoneID <- extract( x=catB, y=zones[[zz]], cellnumbers=TRUE)
  propsOfbathy <- table( catB@data@values[zoneID[[1]][,"cell"]])
  propsOfbathy <- propsOfbathy / sum( propsOfbathy)
  if( zz == "NPZ")
  tmp <- Bathy.targetProps.npz / propsOfbathy
  if( zz == "HPZ")
  tmp <- Bathy.targetProps.hpz / propsOfbathy
  if( zz == "GBS")
  tmp <- Bathy.targetProps / propsOfbathy #the desired inclusion probs (unstandardised)
  for( ii in 1:length( propsOfbathy)){
    inclProbs[zoneID[[1]][,"cell"]][zoneID[[1]][,"value"]==ii] <- tmp[ii]
  }
  inclProbs[zoneID[[1]][,"cell"]][is.na( inclProbs[zoneID[[1]][,"cell"]])] <- 0
  inclProbs[zoneID[[1]][,"cell"]] <- inclProbs[zoneID[[1]][,"cell"]] / sum( inclProbs[zoneID[[1]][,"cell"]])
}
inclProbs@data@values[inclProbs@data@values %in% c(0,1,2,3,4,5,6,7,8)] <- NA  #cheats way to crop
plot( inclProbs)
ip <- inclProbs

#standardising so that the zone totals are correct according to straw.props | straw.nums
hpz <- extract( x=catB, y=zones$HPZ, cellnumbers=TRUE)
npz <- extract( x=catB, y=zones$NPZ, cellnumbers=TRUE)
gbs <- extract( x=catB, y=zones$GBS, cellnumbers=TRUE)

#HPZZone <- extract( x=catB, y=zones$HPZ, cellnumbers=TRUE)
#NPZZone <- extract( x=catB, y=zones$NPZ, cellnumbers=TRUE)

inclProbs@data@values[hpz[[1]][,'cell']] <- inclProbs@data@values[hpz[[1]][,'cell']] * straw.props["HPZ"]
inclProbs@data@values[npz[[1]][,'cell']] <- inclProbs@data@values[npz[[1]][,'cell']] * straw.props["NPZ"]
inclProbs@data@values[gbs[[1]][,'cell']] <- inclProbs@data@values[gbs[[1]][,'cell']] * straw.props["GBS"]

#inclProbs@data@values[HPZZone[[1]][,'cell']] <- inclProbs@data@values[HPZZone[[1]][,'cell']] * straw.props["HPZ"]
#inclProbs@data@values[NPZZone[[1]][,'cell']] <- inclProbs@data@values[NPZZone[[1]][,'cell']] * straw.props["NPZ"]

plot(inclProbs)


writeRaster(inclProbs,paste(d.dir, 'inclProbs_forBOSS_5.tif', sep='/'), overwrite=TRUE)


# to check if sum of cells (probs) = 1
sumr <- cellStats(inclProbs, 'sum')
sumr
