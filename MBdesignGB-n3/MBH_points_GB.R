#############################
### MBH with Points #########

### Libraries ###
library( rgdal)
library( rgeos)
library( sp)
library( raster)
library( MBHdesign)
library( fields)
library( graphics)
library( grDevices) # for colours
library( plot.matrix) # to plot matrices using plot function
library( measurements)

rm( list=ls())

w.dir <- "~/MBHdesignGB"
s.dir <- "~/MBHdesignGB/SpatialData" # spatial data folder
o.dir <- "~/MBHdesignGB/outputs"
d.dir <- "~/MBHdesignGB/data"
p.dir <- "~/MBHdesignGB/plots"


## load bathymetry ----
bathy <- raster(paste(s.dir, "GB_CMR_bathy_utm.tif", sep='/'))
plot(bathy)
proj4string(bathy) # "+proj=utm +zone=50 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
bathy # 36395 cells

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

### 2. Sample acording to depth w MBH package ####
#### Now try the spatially balanced design with preference for High TPI regions #####


# par is a graphics functon used to query graphical parameters
# mfrow = A vector of the form c(nr, nc). Subsequent figures will be drawn in an nr-by-nc array on the 
# device by columns (mfcol), or rows (mfrow), respectively.
# mar: A numerical vector of the form c(bottom, left, top, right) which gives the number of lines 
# of margin to be specified on the four sides of the plot. The default is c(5, 4, 4, 2) + 0.1
par ( mfrow= c ( 1 , 3 ), mar= rep ( 4 , 4 )) 

# set n 
n <- 10

# The number of 'depth bins' to spread sampling effort over.
nbins <- 4

# force the breaks so R doesn't use 'pretty'

breaks <- seq ( from= min ( bathydf$depth, na.rm= TRUE ),
                to= max ( bathydf$depth, na.rm= TRUE ), length= nbins +1 ) # -6.3675003 -0.8616668  4.6441666 10.1500001

# chech the values above with raster
minValue(bathy$GB_CMR_bathy_utm) # -46.90025
maxValue(bathy$GB_CMR_bathy_utm) # -8.540641

# Find sensible tpi bins using pre-packaged code
tmpHist <- hist ( bathydf$depth, breaks= breaks, plot= FALSE )

#  check how it looks to see les frequent values 
# to see the plot
tmpHist <- hist ( bathydf$depth, breaks= breaks, plot= T, freq = F )


# change breaks if needed - least interest areas should have more counts

# When reset the breaks run the tmpHist again
breaks <- c(-46.90025, -30, -20, -8.540641) #
#breaks <- c(-1.1969223, -0.8,  0.5258718,  1.9372688)
# check breaks
tmpHist <- hist ( bathydf$depth, breaks= breaks, freq= F, plot= T )

# Find the inclusion probability for each 'stratum' (for earch 'bin')
tmpHist$inclProbs <- (n/(nbins)) / tmpHist$counts # 1.136364e-01 4.935784e-06 1.138408e-04 3.921569e-01 2.000000e+01

# Matching up locations to probabilties - in data frame
str(bathydf)
tmpHist$ID <- findInterval ( bathydf$depth, tmpHist$breaks) # breaks coded as 1, 2, 3 and so on depending on how many bins
tmpHist$ID[20000] # 3
# not sure why the NAs, but with depth it worked fine
length(tmpHist$ID)

# see hist - quicker way to observe freq of each bin
hist(tmpHist$ID)

head(bathydf)
# A container for the design
# create data frame with each location and its inlcusion probability ####
design <- data.frame ( siteID= 1 : nrow ( bathydf),
                       Easting= bathydf$Easting, Northing= bathydf$Northing,
                       depth= bathydf$depth, inclProb= tmpHist$inclProbs[tmpHist$ID])
str(design)
head(design)


### test ####
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
### test ####


# order design ascending
#design <- design[ order (  design$Easting, design$Northing),]
#head(design)

# remove NA's from design
#design <- design[complete.cases(design),]

# Plot the tpis and the inclusion probabilties 
# toplot this, I need design df
# THIS PLOT TAKES A LONG TIME - find a better way to plot this
#with ( design, plot ( tpi, inclProb, main= "Inclusion Probabilities" ,
#ylab= "Inclusion Probabilities" , xlab= "tpi" , pch= 20 , cex= 1.4 ))

#plot ( design$Depth, design$inclProb, main= "Inclusion Probabilities" ,
#ylab= "Inclusion Probabilities" , xlab= "tpi" , pch= 20 , cex= 1.4 )

# Plot the inclusion probabilities as map
# first, create matrix out of design df
str(design)
length(design$Easting)
# make matrix with inclProbs from design
m <- matrix ( design$inclProb, nrow= length ( uniqueEast), byrow= F)
str(m)
head(m)
m[1]

# then plot
with ( design,
       image.plot ( uniqueEast, uniqueNorth, m,
                    xlab= "" , ylab= "" , main= "Inclusion Probability" , asp= 1 ,
                    ylim= NLims, xlim= ELims))

# Take the Sample using the inclusion probabilities ####
#design$inclProb[4174928] <- 4.935784e-06 ### give this NA an inclusion value

str(design)

##### replace Nas of Inclusion probabilities for zeroes ####
names(design)
head(design)
design$inclProb[is.na(design$inclProb)] <- 0
head(design)
class(design)
any(is.na(design$inclProb))


# turn design df into matrix
designMat <- cbind(design$Easting,design$Northing, design$inclProb)
head(designMat)
str(designMat)
class(designMat)
colnames(designMat) <- c("Easting", "Northing", "inclProbs")
#Fix up ordering issue: revert ordering of columns: first column last so: inclProbs, Northing, Easting
designMat <- designMat[, rev ( 1 : ncol (designMat))]
###############

## Make data frame out of designMat to make raster ###
designDF <- as.data.frame(designMat)
head(designDF)
designDF <- designDF[ order ( designDF$Northing, designDF$Easting),] # order ascending first by northing and then by easting


# turn data frame into raster ####
# df into spatial points 
coordinates(designDF) <- ~ Easting + Northing
# coerce to SpatialPixelsDataFrame
gridded(designDF) <- TRUE
# coerce to raster
designr <- raster(designDF)
designr
plot(designr)


#second Mat - for plotting
designMat2 <- raster::as.matrix(designr, mode ='any')
dim(designMat2)
str(designMat2)
# transpose matrix
designMat3 <- t(designMat2)
dim(designMat3)
str(designMat3)
designMat3 <- designMat3[, rev ( 1 : ncol (designMat3))]

# to see NAs
designr2 <- designr
plot(designr2)
designr2[designr2$inclProbs == 0] <- 20 # give zeroes a different value
plot(designr2)


### make a new data frame out of this raster and Matrix
designdf <- as.data.frame (
  cbind ( coordinates ( designr), as.numeric ( designMat3))) 
colnames ( designdf) <- c ( "Easting" , "Northing" , "inclProbs" ) 
head(designdf)
designdf <- designdf[ order ( designdf$Northing, designdf$Easting),] # order ascending first by northing and then by easting


# Sample with 'quasiSamp' from MBH package #### this takes some time
samp <- quasiSamp ( n = n, dimension= 2 ,
                    potential.sites = coordinates(designr),
                    inclusion.probs= designdf$inclProbs , nSampsToConsider= 100 *n) # inclProb that are not NA!
head(samp)
class(samp)

## Check bathy of the sampling locations
sampsp <- samp
coordinates(sampsp) <- ~x+y
plot(sampsp, add=T)
sampbathy <- raster::extract(bathy, sampsp, df = T)
head(sampbathy)
# check how many samples per "stratum" inclusion probability
sampbathy2 <- cbind(samp, sampbathy)
head(sampbathy2)
str(sampbathy2)
sampbathy2$inclusion.probabilities <- as.factor(sampbathy2$inclusion.probabilities)
summary(sampbathy2)


#Plot the design over bathymetry

# this one is with the last Matrix
#with ( design, image.plot ( uniqueEast, uniqueNorth, designMat3,
#                           xlab= "" , ylab= "" , main= "Spatially-Balanced Sample" , asp= 1 ,
#                          ylim= NLims, xlim= ELims,
#                         col= rev ( tim.colors ())))

image.plot ( uniqueEast, uniqueNorth, bathym2,
             xlab= "Easting" , ylab= "Northing" , main= " " ,
             legend.lab= "Depth" , asp=1 , ylim= NLims, xlim= ELims,
             col=  ( tim.colors ())) 

points ( samp[, c ( "x" , "y" )], pch=20 , cex=2, col = "black" )
