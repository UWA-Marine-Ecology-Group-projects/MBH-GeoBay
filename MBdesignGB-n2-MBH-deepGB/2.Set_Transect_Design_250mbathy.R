library( parallel)
library( class)
library( fields)
library( pdist)
library( raster)
library( rgdal)
library( MBHdesign)

## load bathymetry ----
bathy <- raster(paste(s.dir, "GB250mBathyDeep.tif", sep='/'))
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

#### Set up transect design ####
incprobm <- as.matrix(inclProbs) # need to rearrange the order of this for better plotting
n.x <- nrow( incprobm)
n.y <- ncol( incprobm)
image.plot( x=1:n.x, y=1:n.y, z=incprobm, main="Geographe Bay bathymetry (m)", asp=1)

# transpose 
incprobm2 <- t(incprobm)
n.x <- nrow( incprobm2)
n.y <- ncol( incprobm2)
image.plot( x=1:n.x, y=1:n.y, z=incprobm2, main="Geographe Bay bathymetry (m)", asp=1)

#Fix up ordering issue
incprobm2 <- incprobm2[, rev ( 1 : ncol (incprobm2))]
n.x <- nrow( incprobm2)
n.y <- ncol( incprobm2)
image.plot( x=1:n.x, y=1:n.y, z=incprobm2, main="Geographe Bay bathymetry (m)", asp=1)

cellStats(inclProbs, sum)

#format for MBHdesign functions
incprobdf <- expand.grid( x=1:n.x, y=1:n.y)
incprobdf$IncProbs <- as.vector( incprobm2)
head(incprobdf)


#gb.control is a list that contains transect details
gb.control <- list(
  #the type of transect
  transect.pattern="line",
  #the length of transect
  line.length=250,
  #the number of points that define the transect
  transect.nPts=25,
  #the number of putative directions that a transect can take
  nRotate=11, 
  mc.cores = 25,
  transAdjust = FALSE,
  edge.max.iter = 5,
  nSampsToConsider = 10000
)

## Create constrains, but from previous trials, there are no constrains in GB ####
transAdjust = FALSE### to avoid transects that descend or ascend in bathymetry we use ----
gb.constrains <- findDescendingTrans(
  potential.sites = pot.sites[,c("Easting","Northing")], bathy = pot.sites$depth,
  in.area = rep( TRUE, nrow( pot.sites)), control = gb.control
)

print( dim( gb.constrains)) # [1] 62626    11

#The contents describe how the transect lays over the landscape
#So, there are 15592 putative transects that ascend and descend
# (and can't be used in the sample)
table( as.vector(gb.constrains)) #   563332           121194              410             3950


#convert to TRUE/FALSE
#Note that the final possible transect type ('descendAndNA') is
# not present in these data
#If present, we would have to decide to sample these or not

gb.constraints.bool <- matrix( FALSE, nrow=nrow( gb.constrains),
                               ncol=ncol( gb.constrains))
gb.constraints.bool[gb.constrains %in% c("descend")] <- TRUE

test <- cbind(pot.sites[,c("Easting","Northing")],gb.constraints.bool)
class(test)
coordinates(test) <- ~ Easting + Northing
# coerce to SpatialPixelsDataFrame
gridded(test) <- TRUE
# coerce to raster
constraintsraster <- raster(test)
plot(constraintsraster)

#Let's get a visual to see what has just been done.### this is not working now ----
tmpMat <- matrix( apply( gb.constraints.bool, 1, sum), nrow=nrow( gb.constrains), ncol=ncol( gb.constrains))
tmpMat <- matrix( apply( gb.constraints.bool, 1, sum), nrow=n.y, ncol=n.x)
image.plot( uniqueEast, y=uniqueNorth, z=tmpMat,
            main="Number of Transects",
            sub="Transects centered at cell (max 11)", asp=1)
