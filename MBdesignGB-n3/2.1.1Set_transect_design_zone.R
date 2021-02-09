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
  line.length=1,
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

print( dim( gb.constrains)) # [1] 36395    11

#The contents describe how the transect lays over the landscape
#So, there are 15592 putative transects that ascend and descend
# (and can't be used in the sample)
table( as.vector(gb.constrains)) #  230868    169477 

#convert to TRUE/FALSE
#Note that the final possible transect type ('descendAndNA') is
# not present in these data
#If present, we would have to decide to sample these or not

gb.constraints.bool <- matrix( FALSE, nrow=nrow( gb.constrains),
                               ncol=ncol( gb.constrains))
gb.constraints.bool[gb.constrains %in% c("descend")] <- TRUE

#Let's get a visual to see what has just been done.### this is not working now ----
tmpMat <- matrix( apply( gb.constraints.bool, 1, sum), nrow=nrow( gb.constrains), ncol=ncol( gb.constrains))
tmpMat <- matrix( apply( gb.constraints.bool, 1, sum), nrow=n.x, ncol=n.y)
image.plot( x=1:n.x, y=1:n.y, z=tmpMat,
            main="Number of Transects",
            sub="Transects centered at cell (max 11)", asp=1)


