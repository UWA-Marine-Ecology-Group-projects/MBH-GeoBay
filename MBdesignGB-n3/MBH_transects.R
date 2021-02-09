###########################################

####  MBH transect design for GB #######
## Tnis works with inclusion probabilities 

rm( list=ls())


library( MBHdesign)
library( parallel)
library( class)
library( fields)
library( pdist)
library( raster)
library( rgdal)
library( sp)

w.dir <- "~/MBHdesignGB"
s.dir <- "~/MBHdesignGB/SpatialData" # spatial data folder
o.dir <- "~/MBHdesignGB/outputs"
o2.dir <- "~/MBHdesignGB/outputs/newdesign/"
d.dir <- "~/MBHdesignGB/data"
p.dir <- "~/MBHdesignGB/plots"

## Load GB polygon
gb <- readOGR(paste(s.dir, "GeoBay_CMR_UTM.shp", sep='/'))

## load bathymetry ----
bathy <- raster(paste(s.dir, "GBmultib_lidarUTM_CMR.tif", sep='/'))
plot(bathy)
proj4string(bathy) # "+proj=utm +zone=50 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
bathy # 36395 cells

N <- 36395 # number of cells in domain
#set.seed( 717)
n = 25

bathym <- as.matrix(bathy) # need to rearrange the order of this for better plotting

# conver to matrix for ease plotting
bathym <- raster::as.matrix(bathy)
bathym
str(bathym) # 1:145, 1:251
dim(bathym) # 145 251
bathym[70,75] # -33

# transpose the axis of the matrix so x first and  y later
bathym2 <- t(bathym)
str(bathym2) # [1:251, 1:145]
bathym2[75,70] # -33

# make data frame
bathydf <- as.data.frame (
  cbind (coordinates (bathy), as.numeric (bathym2)))
colnames(bathydf) <- c("Easting", "Northing", "depth")
head(bathydf)
bathydf <- bathydf[ order(bathydf$Northing, bathydf$Easting),] # order ascending first by northing and then by easting
hist(bathydf$depth)

## Setting up plotting for now and later ####
uniqueEast <- base::unique ( bathydf$Easting) # duplicate rows removed
uniqueNorth <- base::unique ( bathydf$Northing)
ELims <- range ( na.exclude ( bathydf)$Easting)
NLims <- range ( na.exclude ( bathydf)$Northing)

str(uniqueEast) ## all the x coordinates
class(uniqueEast)
str(uniqueNorth) ## all the y coordinates
str(bathym2) # the dimensions of the matrix neet to be transposed for the plot
class(bathym2)

#Fix up ordering issue
bathym2 <- bathym2[, rev ( 1 : ncol (bathym2))] # this is needed so the map looks the right way- because of the trasnposing of matrix

## plot it to see what we are dealing with ####
## these kind of plots are kind of slow ###
image.plot ( uniqueEast, uniqueNorth, bathym2,
             xlab= "Easting" , ylab= "Northing" , main= "Geographe Bay CMR" ,
             legend.lab= "Depth" , asp=1 , ylim= NLims, xlim= ELims,
             col=  ( tim.colors ()))     


pot.sites <- bathydf
head(pot.sites)

n.x <- nrow( bathym)
n.y <- ncol( bathym)

#image.plot( x=1:n.x, y=1:n.y, z=bathym, main="Geographe Bay bathymetry (m)", asp=1)
#format for MBHdesign functions
#pot.sites <- expand.grid( x=1:n.x, y=1:n.y)
#pot.sites$depth <- as.vector( bathym)
#head(pot.sites)
#class(pot.sites)



#### Inclusion probabilities ----

# set n 
n <- 20

# The number of 'depth bins' to spread sampling effort over.
nbins <- 4

# force the breaks so R doesn't use 'pretty'

breaks <- seq ( from= min ( bathydf$depth, na.rm= TRUE ),
                to= max ( bathydf$depth, na.rm= TRUE ), length= nbins +1 ) # -46.900249 -34.113713 -21.327177  -8.540641



# chech the values above with raster
minValue(bathy$GB_CMR_bathy_utm) # -46.90025
maxValue(bathy$GB_CMR_bathy_utm) # -8.540641

# Find sensible tpi bins using pre-packaged code
tmpHist <- hist ( bathydf$depth, breaks= breaks, plot= FALSE )

#  check how it looks to see les frequent values 
# to see the plot
tmpHist <- hist ( bathydf$depth, breaks= breaks, plot= T, freq = F )

# Find the inclusion probability for each 'stratum' (for earch 'bin')
tmpHist$inclProbs <- (n/(nbins)) / tmpHist$counts # 1.136364e-01 4.935784e-06 1.138408e-04 3.921569e-01 2.000000e+01

###
#### change breaks if needed - least interest areas should have more counts ----

# When reset the breaks run the tmpHist again

breaks <- c(-46.90025, -30, -20,  -8.540641)

# check breaks
tmpHist <- hist ( bathydf$depth, breaks= breaks, freq= F, plot= T )

# Find the inclusion probability for each 'stratum' (for earch 'bin')
tmpHist$inclProbs <- (n/(nbins)) / tmpHist$counts


# Matching up locations to probabilties - in data frame
str(bathydf)
tmpHist$ID <- findInterval ( bathydf$depth, tmpHist$breaks) # breaks coded as 1, 2, 3 and so on depending on how many bins
tmpHist$ID[20000] # 3
# not sure why the NAs, but with depth it worked fine
length(tmpHist$ID)

# see hist - quicker way to observe freq of each bin
hist(tmpHist$ID)
hist(tmpHist$inclProbs)

head(bathydf)
# A container for the design
# create data frame with each location and its inlcusion probability ####
design <- data.frame ( siteID= 1 : nrow ( bathydf),
                       Easting= bathydf$Easting, Northing= bathydf$Northing,
                       depth= bathydf$depth, inclProb= tmpHist$inclProbs[tmpHist$ID])
str(design)
head(design)

incprob <- design[,c(2,3,5)]
head(incprob)

# make df a raster --
coordinates(incprob) <- ~Easting+Northing
gridded(incprob) <- TRUE
rasterIP <- raster(incprob)
plot(rasterIP)
cellStats(rasterIP, sum)

## for function --

incprobm <- as.matrix(rasterIP) # need to rearrange the order of this for better plotting
n.x <- nrow( incprobm)
n.y <- ncol( incprobm)
image.plot( x=1:n.x, y=1:n.y, z=incprobm, main="Geographe Bay bathymetry (m)", asp=1)
#format for MBHdesign functions
incprobdf <- expand.grid( x=1:n.x, y=1:n.y)
incprobdf$IncProbs <- as.vector( incprobm)
head(incprobdf)



#gb.control is a list that contains transect details
gb.control <- list(
  #the type of transect
  transect.pattern="line",
  #the length of transect
  line.length=1000,
  #the number of points that define the transect
  transect.nPts=20,
  #the number of putative directions that a transect can take
  nRotate=11, 
  mc.cores = 25,
  transAdjust = FALSE,
  edge.max.iter = 10,
  nSampsToConsider = 10000
)


transAdjust = FALSE### to avoid transects that descend or ascend in bathymetry we use ----
gb.constrains <- findDescendingTrans(
  potential.sites = pot.sites[,c("Easting","Northing")], bathy = pot.sites$depth,
  in.area = rep( TRUE, nrow( pot.sites)), control = gb.control
)
# this is a matrix with nrow given by the number of sites and ncol by
# the number of rotations around each site

print( dim( gb.constrains)) # [1] 36395    11

#The contents describe how the transect lays over the landscape
#So, there are 15592 putative transects that ascend and descend
# (and can't be used in the sample)
table( as.vector(gb.constrains)) #  224105            63664            10978           101598  

#convert to TRUE/FALSE
#Note that the final possible transect type ('descendAndNA') is
# not present in these data
#If present, we would have to decide to sample these or not

gb.constraints.bool <- matrix( FALSE, nrow=nrow( gb.constrains),
                                ncol=ncol( gb.constrains))
gb.constraints.bool[gb.constrains %in% c("descend")] <- TRUE

#Let's get a visual to see what has just been done.### this is not working now ----
tmpMat <- matrix( apply( gb.constraints.bool, 1, sum), nrow=nrow( gb.constrains), ncol=ncol( gb.constrains))
tmpMat <- matrix( apply( gb.constraints.bool, 1, sum), nrow=n.y, ncol=n.x)
image.plot( uniqueEast, uniqueNorth, z=tmpMat,
            main="Number of Transects",
            sub="Transects centered at cell (max 11)", asp=1)
#There aren't any transects that are centred on ridges or depressions


#take the sample
GBSamp <- transectSamp( n=n, potential.sites=pot.sites[,c("Easting","Northing")],
                        inclusion.probs= incprobdf[,3],
                         control=gb.control,
                        constrainedSet=gb.constraints.bool
                         )

#visualise the sample

image.plot ( uniqueEast, uniqueNorth, bathym2,
             xlab= "Easting" , ylab= "Northing" , main= "Geographe Bay CMR" ,
             legend.lab= "Depth" , asp=1 , ylim= NLims, xlim= ELims,
             col=  ( tim.colors ()))
points( GBSamp$points[,c("Easting","Northing")], pch=20)

#write csv
write.csv( GBSamp$transect, paste(o2.dir, "1kmSpatBal20Transects.csv", sep='/'), row.names=FALSE)
write.csv(GBSamp$points, paste(o2.dir, "1kmSpatBal20Transectpoints.csv", sep='/'), row.names=FALSE)
#tidy

Tpoints <- GBSamp$points
coordinates(Tpoints) <- ~Easting+Northing

pdf( "~/MBHdesignGB/outputs/newdesign/1kmSpatBal20Transects.pdf", height=7, width=8)
plot( bathy, main="Spatially balanced transects 1km (n=20)")
points(Tpoints, add=T, pch=20, col = 'black', add = TRUE)
plot( gb, add=TRUE)
dev.off()

writeOGR(Tpoints, dsn="~/MBHdesignGB/outputs/newdesign", layer=paste( "1kmSpatBal20Transects", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)

rm( list=ls())





###################
####   PLOT   ####
##################

# Plot bathy
proj4string(bathy) # "+proj=utm +zone=50 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(bathy)

# load and plot CMR zones 
cmr <- readOGR(paste(s.dir, "GeoBay_CMR_UTM.shp", sep='/'))
plot(cmr, add=T)
proj4string(cmr)

# Plot transect locations
head(GBSamp$transect)
samp <- GBSamp$transect
head(samp)
coordinates(samp) <- ~Easting+Northing
points(samp, add=T)
