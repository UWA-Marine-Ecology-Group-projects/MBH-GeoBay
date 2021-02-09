### Create inclusion probabilities ####

## for multibeam and lidar #######

library( MBHdesign)
library( parallel)
library( class)
library( fields)
library( pdist)
library( raster)
library( rgdal)
library( sp)

#########################
#read in data

setwd("~/MBHdesignGB/")

GBDat <- readRDS( "GBData_forDesign3.RDS")
gb_rasters <- readRDS("GBRasters_forDesign3.RDS")
zones <- readRDS( "GBZones_forDesign3.RDS")

####################################
####  Straw man for numbers of samples in each region
####################################

straw.nums <- c( 16, 16, 6)  # numbers of drops
straw.props <- straw.nums / sum( straw.nums)
names( straw.nums) <- names( straw.props) <- c( "MUZ", "SPZ", "HPZ")
saveRDS( straw.nums, file="StrawmanNumbers_Zones.RDS")

setwd("~/MBHdesignGB/")


###################################
####  Hand-picking Bathy cut points
####  And their numbers of drops
###################################

#Bathy.quant <- c(0,0.8,0.9,0.925,0.95,0.975,1)
Bathy.quant <- c(0,0.25,0.5,0.8,1)
Bathy.cuts <- quantile( gb_rasters$bathy, Bathy.quant)#c( -Inf,0.02,0.04,0.08,0.16,Inf)
#trying to make it so there is no hand-picking (except for the hand-picked function)
tmp <- cumsum( Bathy.quant)
Bathy.targetNums <- rep( floor( 18/8), 4)#floor( ( tmp / sum( tmp))[-1] * 200)#rep( 40, 5)#c( 20,20,30,65,65)
Bathy.targetProps <-  Bathy.targetNums / sum( Bathy.targetNums)

########################
#depth limiting and some elementary cleaning

#minDepthInAMP <- max( NingalooDat$BATHY, na.rm=TRUE)

#NingalooDat[ is.na( NingalooDat$BATHY) | NingalooDat$BATHY < -195 & NingalooDat$BATHY > minDepthInAMP, c("BATHY","TPI_GF")] <- NA

#GBDat[ !is.na( GBDat$BATHY) & GBDat$BATHY < -50 , "BATHY"] <- NA
#GBDat[ !is.na( GBDat$BATHY) & GBDat$BATHY > 0 , "BATHY"] <- NA


pdf( "Bathy_distribution_zones_MultibLidar.pdf")
{
  
  par( mfrow=c(2,2))
  for( ii in c("MUZ", "SPZ", "HPZ", "NPZ"))
    hist( GBDat[GBDat[,ii]==1,"BATHY"], nclass=50, main=ii)
  
  par( mfrow=c(1,1))
  kount <- 1
  plot( ecdf( GBDat[GBDat[,"MUZ"]==1,"BATHY"]), col=2, main="Bathy in Zones", xlim=c(-0.4,0.4))
  for( ii in c("MUZ", "SPZ", "HPZ", "NPZ")){
    plot( ecdf( GBDat[GBDat[,ii]==1,"BATHY"]), col=kount, add=TRUE)
    kount <- kount + 1
  }
  legend( "bottomright", legend=c("MUZ", "SPZ", "HPZ", "NPZ"), lty=1, 
          col=c(1:4), lwd=c(1,1,1,1), bty='n')
  #so surprisingly, the southern control zone has more flat area than the IUCN2 area
  #there may be some issues with the north control having quite a few more 'bumps' than the southern control and IUCN2 area.  Bias.
  dev.off()
}
#########################
#proportion of potential sites in each zone

GBDat_small <- GBDat[!is.na( GBDat$BATHY),]
tmp <- colSums( GBDat_small[,c("MUZ", "SPZ", "HPZ")], na.rm=TRUE) 
#tmp[2] <- tmp[2] - tmp[1] # so similar amount of sites in SPZ and MUZ
props <- tmp / nrow( GBDat_small)
props <- props / sum( props) # 1 UP TO HERE

###################################
####  TPI to get cut points
###################################

catB <- cut( gb_rasters$bathy, breaks=Bathy.cuts, na.rm=TRUE)

plot( zones$MUZ); plot( catB, add=TRUE); plot( zones$MUZ, add=TRUE)

writeRaster( catB, file='Bathy_cuts_MultibLidar.tif', overwrite=TRUE)

plot(catB)

##################################
####  Within each zone (incl probs)
####  Weight according to straw.props
##################################

Bathy.targetProps <- c(0.25,0.25,0.25,0.25)
Bathy.targetProps2 <- c(0.5,0.5)

inclProbs <- catB
for( zz in c("HPZ", "SPZ", "MUZ")){
  print( zz)
  #if( zz == "MUZ")
  #zoneID <- extract( x=catB, y=zones$MUZ, cellnumbers=TRUE)
  #zoneID <- extract( x=catB, y=zones$MUZ-zones$NPZ, cellnumbers=TRUE)
  #else
  zoneID <- extract( x=catB, y=zones[[zz]], cellnumbers=TRUE)
  propsOfTPI <- table( catB@data@values[zoneID[[1]][,"cell"]])
  propsOfTPI <- propsOfTPI / sum( propsOfTPI)
  if(length(propsOfTPI) == 4)
    tmp <- Bathy.targetProps / propsOfTPI #the desired inclusion probs (unstandardised)
  else
    tmp <- Bathy.targetProps2 / propsOfTPI #the desired inclusion probs (unstandardised)
  for( ii in 1:length( propsOfTPI)){
    inclProbs[zoneID[[1]][,"cell"]][zoneID[[1]][,"value"]==ii] <- tmp[ii]
  }
  inclProbs[zoneID[[1]][,"cell"]][is.na( inclProbs[zoneID[[1]][,"cell"]])] <- 0
  inclProbs[zoneID[[1]][,"cell"]] <- inclProbs[zoneID[[1]][,"cell"]] / sum( inclProbs[zoneID[[1]][,"cell"]])
}
inclProbs@data@values[inclProbs@data@values %in% c(0,1,2,3,4,5,6,7,8)] <- NA  #cheats way to crop
plot( inclProbs)

#standardising so that the zone totals are correct according to straw.props | straw.nums
MUZzone <- extract( x=catB, y=zones$MUZ, cellnumbers=TRUE)
SPZzone <- extract( x=catB, y=zones$SPZ, cellnumbers=TRUE)
HPZZone <- extract( x=catB, y=zones$HPZ, cellnumbers=TRUE)
#NPZZone <- extract( x=catB, y=zones$NPZ, cellnumbers=TRUE)

inclProbs@data@values[MUZzone[[1]][,'cell']] <- inclProbs@data@values[MUZzone[[1]][,'cell']] * straw.props["MUZ"]
inclProbs@data@values[SPZzone[[1]][,'cell']] <- inclProbs@data@values[SPZzone[[1]][,'cell']] * straw.props["SPZ"]
inclProbs@data@values[HPZZone[[1]][,'cell']] <- inclProbs@data@values[HPZZone[[1]][,'cell']] * straw.props["HPZ"]
#inclProbs@data@values[NPZZone[[1]][,'cell']] <- inclProbs@data@values[NPZZone[[1]][,'cell']] * straw.props["NPZ"]

plot(inclProbs)

writeRaster( inclProbs, file='inclProbs_design4.tif', overwrite=TRUE)

rm( list=lso()$OTHERS)
rm( AMPzone, catTPI, ii, IUCNzone, kount, minDepthInAMP, NingalooDat_small, NthZone, props, propsOfTPI, SthZone, 
    tmp, zoneID, zz, TPI.cuts, TPI.quant, TPI.targetNums, TPI.targetProps)

