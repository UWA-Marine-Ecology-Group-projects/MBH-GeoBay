###################################################################################################
#####   Simulation to investigate properties of transect samples
#####	Inclusion probability map from the volcano data set
#####   Written by Scott Aug-18
###################################################################################################

####################################################################
####	Run on R-3.6.1 during July 2019
####	We note that the random number generator was updated for 
####	R-3.6.x and so results may differ if script is run on a previous
####	version of R.
####################################################################

#########################################################
####  Clean and set up libraries first
#########################################################

rm( list=ls())
#.libPaths( c(.libPaths(), "~/lib/R/library"))

library( MBHdesign)
library( parallel)
library( class)
library( fields)
library( pdist)

#########################################################
####  Some functions for choosing transects naively
#########################################################

#Choose a transect direction given its centre location (naive)
"niave_chooseDirection" <- function( centreLoc, designParams, transectParams){
  directSet <- 1:designParams$N$rotate
  direct <- sample( x=directSet, size=1)
  flag <- TRUE
  while( flag){
    tran <- matrix( rep( as.numeric( centreLoc), each=designParams$N$transect.nPts), ncol=2) + transectParams$potential.transects[[direct]]
    inside <- mgcv:::in.out(designParams$study.area, tran)
    if( all( inside))
      flag <- FALSE
    else{
      directSet <- setdiff(directSet, direct)
      if( length( directSet)==0)
        return( NA)
      direct <- sample( x=directSet, size=1)
    }
  }
  return( direct)
}

#Design using Naive methods (SRS or BAS based on cell probabilities)
"sampleTransect_withCellProbs" <- function( x, n, inclusion.probs=NULL, potential.sites=NULL, 
                                            study.area=NULL, control=NULL){
  #argument x is ignored.
  suppressMessages( control <- MBHdesign:::set.transect.control( control))
  #set the parameters of the design (survey area, potential sites, inclusion probs, all the different n's, ...)
  suppressMessages( designParams <- MBHdesign:::setDesignParams_transect( study.area=study.area, potential.sites=potential.sites, inclusion.probs=inclusion.probs, control=control))
  #set parameters of the individual transects (pattern, direction, point representation, ...)
  suppressMessages( transectParams <- MBHdesign:::setTransectParams2( designParams=designParams, control=control))
  
  oversampMult <- 10
  totSamp <- oversampMult * n
  if( control$spat.random.type == "pseudo"){
    centreLocs <- sample( x=1:designParams$N$Tot.xy, size=totSamp, replace=FALSE, prob=designParams$inclusion.probs)
    centreLocs <- designParams$potential.sites[centreLocs,]
  }
  if( control$spat.random.type == "quasi")
    suppressMessages( centreLocs <- MBHdesign:::quasiSamp( n=totSamp, dimension=2, study.area=designParams$study.area, potential.sites=designParams$potential.sites, 
                                                           inclusion.probs=designParams$inclusion.probs, randStartType=control$randStartType, nSampsToConsider=control$nSampsToConsider))
  
  
  directions <- rep( NA, totSamp)  
  keepgoing <- TRUE
  kount <- 1
  while( keepgoing==TRUE){
    directions[kount] <- niave_chooseDirection( centreLocs[kount,1:2], designParams, transectParams)
    kount <- kount + 1
    if( sum( !is.na( directions)) == n)
      keepgoing <- FALSE
  }
  ids <- which( !is.na( directions))
  funny1 <- function(x){
    matrix( rep( as.numeric( centreLocs[x,1:2]), each=designParams$N$transect.nPts), ncol=2) + transectParams$potential.transects[[directions[x]]]
  }
  pts.list <- lapply( ids, FUN=funny1)
  pts <- do.call( "rbind", pts.list)
  
  return( pts)  
}

##########################################################
####	Some functions for choosing stratified transects
##########################################################

"findStrata_fromInclProbs" <- function( x, nstrata = 5){
  minny <- min( x)
  maxxy <- max( x)
  cumSummy <- rep( NA, length( unique( x)+1))
  kount <- 1
  for( ii in sort( unique( x))){
    cumSummy[kount] <- sum( x[x<=ii])
    kount <- kount + 1
  }
  cumSummy <- cumSummy / max( cumSummy)
  cuts.ids <- rep( NA, 1+nstrata)
  for( ii in 0:(nstrata)){
    cuts.ids[ii+1] <- which.min( abs( cumSummy - (max( cumSummy) - min( cumSummy)) * ii / nstrata))
  }
  
  strata <- findInterval( x, sort( unique( x))[cuts.ids])
  strata[strata==nstrata+1] <- nstrata
  
  #image( array( strata, dim=dim( x)))
  
  return( strata)
}

#Design using Naive stratification methods (stratified on cell probs -- SRS only)
"sampleTransect_strat" <- function( x, n, strata, inclusion.probs=NULL, potential.sites=NULL, 
                                    study.area=NULL, control=NULL){
  #argument x is ignored.
  suppressMessages( control <- MBHdesign:::set.transect.control( control))
  #set the parameters of the design (survey area, potential sites, inclusion probs, all the different n's, ...)
  suppressMessages( designParams <- MBHdesign:::setDesignParams_transect( study.area=study.area, potential.sites=potential.sites, inclusion.probs=inclusion.probs, control=control))
  #  #set parameters of the individual transects (pattern, direction, point representation, ...)
  suppressMessages( transectParams <- MBHdesign:::setTransectParams2( designParams=designParams, control=control))
  
  #set the numbers of transects per strata
  nstrata <- length( unique( strata))
  nperstrata <- rep( n %/% nstrata, nstrata)
  oversampMult <- 10
  #the left over transects
  if( sum( nperstrata) != n){
    nremain <- n %% nstrata
    increaseStrata <- sample( 1:nstrata, size=nremain, replace=FALSE)
    nperstrata[increaseStrata] <- nperstrata[increaseStrata] + 1
  }
  #centre IDs, one for each strata
  centreIDs <- centreLocs <- list()
  for( ii in 1:nstrata){
    centreIDs[[ii]] <- sample( x=(1:designParams$N$Tot.xy)[strata==ii], size=oversampMult * nperstrata[ii], replace=FALSE)
    centreLocs[[ii]] <- designParams$potential.sites[centreIDs[[ii]],]
  }
  #the directions for those locations, if possible.
  directions <- list()
  for( ii in 1:nstrata){
    keepgoing <- TRUE
    kount <- 1
    directions[[ii]] <- rep( NA, oversampMult * nperstrata[ii])
    while( keepgoing==TRUE){
      directions[[ii]][kount] <- niave_chooseDirection( centreLocs[[ii]][kount,1:2], designParams, transectParams)
      kount <- kount + 1
      if( sum( !is.na( directions[[ii]])) == nperstrata[ii])
        keepgoing <- FALSE
    }
  }
  ids <- lapply( directions, function(x) which( !is.na( x)))
  funny1 <- function(x, clocs, dirs){
    matrix( rep( as.numeric( clocs[x,1:2]), each=designParams$N$transect.nPts), ncol=2) + transectParams$potential.transects[[dirs[x]]]
  }
  pts.list <- list()
  for( ii in 1:nstrata){
    pts.list[[ii]] <- lapply( ids[[ii]], FUN=funny1, clocs=centreLocs[[ii]], dirs=directions[[ii]])
  }
  pt.list1 <- lapply( pts.list, function(x) do.call( "rbind", x))
  pts <- do.call( "rbind", pt.list1)
  
  return( pts)  
}

#########################################################
####  A function for repeated sampling of transect-designs
#########################################################

funny <- function(x, control){
  control1 <- control
  control1$return.index <- FALSE
  control1$calcObsProbs <- FALSE
  control$mc.cores <- 1
  #  suppressMessages( tmp <- transectSamp(n=n, study.area=studyArea, potential.sites=X, inclusion.probs=p, control=control1,
  #                                        IDs=generalPattern$index$IDs, adjustedSpecified=generalPattern$index$adjustedSpecified))
  suppressMessages( tmp <- MBHdesign:::transectSamp.internal(n=n, study.area=studyArea, potential.sites=X, inclusion.probs=p, control=control1,
                                                             IDs=generalPattern$index$IDs, adjustedSpecified=generalPattern$index$adjustedSpecified))
  ret <- list( transectStarts = tmp$transect[,c("x","y","direction")], points = tmp$points[,c("x","y")])
  if( x %% 1000 == 0)
    cat( x, " ")
  return( ret)
}

#########################################################
####  A function for calculating transect distances from
###   a matrix with a particular form
#########################################################

"designTranDist" <- function( kk){
  #calculate the average transect distance, and its sd
  distMatMin <- distMatMean <- matrix( NA, ncol=n, nrow=n)
  for( ii in 1:n){
    if( ii < n)
      for (jj in (ii+1):n){
        tran1 <- manyDesigns[[kk]][(ii-1)*tranNp + 1:tranNp,]
        tran2 <- manyDesigns[[kk]][(jj-1)*tranNp + 1:tranNp,]
        distMatMin[ii,jj] <- min( pdist( tran1, tran2)@dist)
        distMatMean[ii,jj] <- mean( pdist( tran1, tran2)@dist)
      }
  }
  ret <- c( min( distMatMin[upper.tri(distMatMin)]), min( distMatMean[upper.tri(distMatMean)]))
  names( ret) <- c("min","mean")
  return( ret)
}

#########################################################
####  Functions to calculate weighted Voronoi variances
#########################################################

distFunny <- function( refPt, pts){
  tmp <- pts - as.matrix( matrix( refPt, ncol=2)[rep( 1,nrow( pts)),])
  tmp1 <- tmp^2
  tmp2 <- rowSums( tmp1)
  tmp3 <- sqrt( tmp2)
  return( tmp3)
}

getVarArea <- function( des, X.grid, inclProbs){
  disty <- apply( des, 1, distFunny, pts=X.grid)
  nearPt <- apply( disty, 1, which.min)
  nearTrans <- (nearPt-1) %/% control$transect.nPts + 1
  inclProbs <- ( inclProbs / sum( inclProbs)) * n
  vi <- tapply( inclProbs, nearTrans, sum)
  
  meany <- mean( vi)
  vary <- sum( (vi-1)^2) / length( vi)
  
  ret <- c(meany,vary)
  names( ret) <- c("mean","var")
  
  return( ret)
}

#########################################################
####  Number of simulations
#########################################################

nSims <- 200000

#########################################################
####  Design parameters
#########################################################

#the number of potential sampling locations
n.x <- nrow( volcano)
n.y <- ncol( volcano)
N <- n.x*n.y

#number of samples
n <- 10
#the grid on unit square
X <- as.matrix( expand.grid( x=10*(1:n.x)-5, y=10*(1:n.y)-5))	#cells are 10mx10m approx
##the inclusion probabilities with gradient according to non-linear function of X[,1] and X[,2]
p <- as.numeric( volcano)

#image( x=10*(1:n.x)-5, y=10*(1:n.y)-5, z=matrix( p, ncol=n.y, nrow=n.x))
#standardise
p <- p / sum( p)
#The study area (unit square)
studyArea <- matrix( c( 0,0, 10*n.x,0, 10*n.x,10*n.y, 0,10*n.y), ncol=2, byrow=TRUE)
#The number of points to represent each transect
tranNp <- 18	#more than one point per cell for horizontal transects
#The number of rotations to consider when randomly sampling
nRot <- 51
#length of linear transect
lineLength <- 175
#Number of cores (max for machine)
mc.cores <- parallel::detectCores()

#number of strata for stratified sampling
nstrata <- 5

#########################################################
####  Containers for storing results
#########################################################

results <- list()
results$spat.bal <- matrix( NA, nrow=6, ncol=5)
results$incl.probs.mse <- matrix( NA, ncol=5, nrow=2)
colnames( results$spat.bal) <- colnames( results$incl.probs.mse) <- c("Transect_pseudo","Transect_quasi","cell_pseudo","cell_quasi","strata")
rownames( results$spat.bal) <- c("mean_var","sd_var", "mean_minTransDist", "sd_minTransDist", "mean_meanTransDist", "sd_meanTransDist")
rownames( results$incl.probs.mse) <- c( "mse","absDev")
cellProps <- list()

#########################################################
####  BAS transects (based on transect probabilities)
####  WITH edge adjustment
#########################################################
{
  rand.type <- "quasi"
  
  #set up design/function parameters
  control <- MBHdesign:::set.transect.control( list( nRotate=nRot, transect.nPts=tranNp, line.length=lineLength,
                                                     edgeAdjust=TRUE, spat.random.type=rand.type, return.index=TRUE))
  
  #get the computationally expensive bits
  generalPattern <- transectSamp(n=n, study.area=studyArea, potential.sites=X, inclusion.probs=p, control=control)
  #get lots of designs
  combinedSamples <- parallel::mclapply( 1:nSims, FUN=funny, mc.cores=mc.cores, control=control)
  #combine them (matrices are apparently quicker)
  manyDesigns <- parallel::mclapply( combinedSamples, FUN=function(x) as.matrix( x$points))
  pointSamples <- matrix( NA, ncol=2, nrow=nSims * n * tranNp)
  for( ii in 1:nSims)
    pointSamples[(ii-1)*n*tranNp + 1:(n*tranNp),] <- manyDesigns[[ii]]
  pointSamples <- as.data.frame( pointSamples)
  
  #find the empirical (frequentist) probability of observing each cell
  pointSamples$cellID <- class::knn1( train=X, test=pointSamples, cl=1:N)
  cellVisits <- table( pointSamples$cellID) #this is already a factor with N levels.
  cellProps$Transect_quasi <- cellVisits / sum( cellVisits)
  
  #calculate the MSE
  results$incl.probs.mse["mse","Transect_quasi"] <- sum((cellProps$Transect_quasi-p)^2)/N
  #and the abs deviation
  results$incl.probs.mse["absDev","Transect_quasi"] <- mean( abs( cellProps$Transect_quasi-p))/N
  
  #calculate the average spat. bal. and its sd
  voronoiVars <- mclapply( manyDesigns, getVarArea, X.grid=X, inclProbs=p, mc.cores=mc.cores)
  voronoiVars <- t( simplify2array( voronoiVars))
  
  results$spat.bal["mean_var","Transect_quasi"] <- mean( voronoiVars[,"var"])
  results$spat.bal["sd_var","Transect_quasi"] <- sd( voronoiVars[,"var"])
  
  tmp <- t( simplify2array( mclapply( 1:nSims, designTranDist, mc.cores=mc.cores)))
  
  results$spat.bal["mean_minTransDist","Transect_quasi"] <- mean( tmp[,'min'])
  results$spat.bal["sd_minTransDist","Transect_quasi"] <- sd( tmp[,'min'])
  results$spat.bal["mean_meanTransDist","Transect_quasi"] <- mean( tmp[,'mean'])
  results$spat.bal["sd_meanTransDist","Transect_quasi"] <- sd( tmp[,'mean'])
}

#########################################################
####  SRS transects (based on transect probabilities)
####  WITH edge adjustment
#########################################################
{
  rand.type <- "pseudo"
  
  #set up design/function parameters
  control <- MBHdesign:::set.transect.control( list( nRotate=nRot, transect.nPts=tranNp, line.length=lineLength,
                                                     edgeAdjust=TRUE, spat.random.type=rand.type, return.index=TRUE))
  
  #get the computationally expensive bits
  generalPattern <- transectSamp(n=n, study.area=studyArea, potential.sites=X, inclusion.probs=p, control=control)
  #get lots of designs
  combinedSamples <- parallel::mclapply( 1:nSims, FUN=funny, mc.cores=mc.cores, control=control)
  #combine them (matrices are apparently quicker)
  manyDesigns <- parallel::mclapply( combinedSamples, FUN=function(x) as.matrix( x$points))
  pointSamples <- matrix( NA, ncol=2, nrow=nSims * n * tranNp)
  for( ii in 1:nSims)
    pointSamples[(ii-1)*n*tranNp + 1:(n*tranNp),] <- manyDesigns[[ii]]
  pointSamples <- as.data.frame( pointSamples)
  
  #find the empirical (frequentist) probability of observing each cell
  pointSamples$cellID <- class::knn1( train=X, test=pointSamples, cl=1:N)
  cellVisits <- table( pointSamples$cellID) #this is already a factor with N levels.
  cellProps$Transect_pseudo <- cellVisits / sum( cellVisits)
  
  #calculate the MSE
  results$incl.probs.mse["mse","Transect_pseudo"] <- sum((cellProps$Transect_pseudo-p)^2)/N
  #and the abs deviation
  results$incl.probs.mse["absDev","Transect_pseudo"] <- mean( abs( cellProps$Transect_pseudo-p))/N
  
  #calculate the average spat. bal. and its sd
  voronoiVars <- mclapply( manyDesigns, getVarArea, X.grid=X, inclProbs=p, mc.cores=mc.cores)
  voronoiVars <- t( simplify2array( voronoiVars))
  
  results$spat.bal["mean_var","Transect_pseudo"] <- mean( voronoiVars[,"var"])
  results$spat.bal["sd_var","Transect_pseudo"] <- sd( voronoiVars[,"var"])
  
  tmp <- t( simplify2array( mclapply( 1:nSims, designTranDist, mc.cores=mc.cores)))
  
  results$spat.bal["mean_minTransDist","Transect_pseudo"] <- mean( tmp[,'min'])
  results$spat.bal["sd_minTransDist","Transect_pseudo"] <- sd( tmp[,'min'])
  results$spat.bal["mean_meanTransDist","Transect_pseudo"] <- mean( tmp[,'mean'])
  results$spat.bal["sd_meanTransDist","Transect_pseudo"] <- sd( tmp[,'mean'])
}

#########################################################
####  SRS transects (based on CELL probabilities)
#########################################################
{
  rand.type <- "pseudo"
  
  #set up design/function parameters
  control <- MBHdesign:::set.transect.control( list( nRotate=nRot, transect.nPts=tranNp, line.length=lineLength,
                                                     edgeAdjust=FALSE, spat.random.type=rand.type, return.index=TRUE, randStartType=2))
  
  #get lots of designs  
  combinedSamples <- parallel::mclapply( 1:nSims, FUN=sampleTransect_withCellProbs, mc.cores=mc.cores, n=n, 
                                         inclusion.probs=p, potential.sites=X, study.area=studyArea, control=control)
  #combine them (matrices are apparently quicker)
  manyDesigns <- parallel::mclapply( combinedSamples, FUN=function(x) as.matrix( x))
  pointSamples <- matrix( NA, ncol=2, nrow=nSims * n * tranNp)
  for( ii in 1:nSims)
    pointSamples[(ii-1)*n*tranNp + 1:(n*tranNp),] <- manyDesigns[[ii]]
  pointSamples <- as.data.frame( pointSamples)
  
  #find the empirical (frequentist) probability of observing each cell
  pointSamples$cellID <- class::knn1( train=X, test=pointSamples, cl=1:N)
  cellVisits <- table( pointSamples$cellID) #this is already a factor with N levels.
  cellProps$cell_pseudo <- cellVisits / sum( cellVisits)
  
  #calculate the MSE
  results$incl.probs.mse["mse","cell_pseudo"] <- sum((cellProps$cell_pseudo-p)^2)/N
  #and the abs deviation
  results$incl.probs.mse["absDev","cell_pseudo"] <- mean( abs( cellProps$cell_pseudo-p))/N
  
  #calculate the average spat. bal. and its sd
  voronoiVars <- mclapply( manyDesigns, getVarArea, X.grid=X, inclProbs=p, mc.cores=mc.cores)
  voronoiVars <- t( simplify2array( voronoiVars))
  
  results$spat.bal["mean_var","cell_pseudo"] <- mean( voronoiVars[,"var"])
  results$spat.bal["sd_var","cell_pseudo"] <- sd( voronoiVars[,"var"])
  
  tmp <- t( simplify2array( mclapply( 1:nSims, designTranDist, mc.cores=mc.cores)))
  
  results$spat.bal["mean_minTransDist","cell_pseudo"] <- mean( tmp[,'min'])
  results$spat.bal["sd_minTransDist","cell_pseudo"] <- sd( tmp[,'min'])
  results$spat.bal["mean_meanTransDist","cell_pseudo"] <- mean( tmp[,'mean'])
  results$spat.bal["sd_meanTransDist","cell_pseudo"] <- sd( tmp[,'mean'])
}

#########################################################
####  BAS transects (based on CELL probabilities)
#########################################################
{
  rand.type <- "quasi"
  #set up design/function parameters
  control <- MBHdesign:::set.transect.control( list( nRotate=nRot, transect.nPts=tranNp, line.length=lineLength,
                                                     edgeAdjust=FALSE, spat.random.type=rand.type, return.index=TRUE, randStartType=2))
  
  #get lots of designs  
  combinedSamples <- parallel::mclapply( 1:nSims, FUN=sampleTransect_withCellProbs, mc.cores=mc.cores, n=n, 
                                         inclusion.probs=p, potential.sites=X, study.area=studyArea, control=control)
  #combine them (matrices are apparently quicker)
  manyDesigns <- parallel::mclapply( combinedSamples, FUN=function(x) as.matrix( x))
  pointSamples <- matrix( NA, ncol=2, nrow=nSims * n * tranNp)
  for( ii in 1:nSims)
    pointSamples[(ii-1)*n*tranNp + 1:(n*tranNp),] <- manyDesigns[[ii]]
  pointSamples <- as.data.frame( pointSamples)
  
  #find the empirical (frequentist) probability of observing each cell
  pointSamples$cellID <- class::knn1( train=X, test=pointSamples, cl=1:N)
  cellVisits <- table( pointSamples$cellID) #this is already a factor with N levels.
  cellProps$cell_quasi <- cellVisits / sum( cellVisits)
  
  #calculate the MSE
  results$incl.probs.mse["mse","cell_quasi"] <- sum((cellProps$cell_quasi-p)^2)/N
  #and the abs deviation
  results$incl.probs.mse["absDev","cell_quasi"] <- mean( abs( cellProps$cell_quasi-p))/N
  
  #calculate the average spat. bal. and its sd
  voronoiVars <- mclapply( manyDesigns, getVarArea, X.grid=X, inclProbs=p, mc.cores=mc.cores)
  voronoiVars <- t( simplify2array( voronoiVars))
  
  results$spat.bal["mean_var","cell_quasi"] <- mean( voronoiVars[,"var"])
  results$spat.bal["sd_var","cell_quasi"] <- sd( voronoiVars[,"var"])
  
  tmp <- t( simplify2array( mclapply( 1:nSims, designTranDist, mc.cores=mc.cores)))
  
  results$spat.bal["mean_minTransDist","cell_quasi"] <- mean( tmp[,'min'])
  results$spat.bal["sd_minTransDist","cell_quasi"] <- sd( tmp[,'min'])
  results$spat.bal["mean_meanTransDist","cell_quasi"] <- mean( tmp[,'mean'])
  results$spat.bal["sd_meanTransDist","cell_quasi"] <- sd( tmp[,'mean'])
}

#########################################################
####  SRS stratified transects
#########################################################
{
  rand.type <- "strat"
  
  #define the strata
  myStrata <- findStrata_fromInclProbs( x=p, nstrata = nstrata)
  
  #set up design/function parameters
  control <- MBHdesign:::set.transect.control( list( nRotate=nRot, transect.nPts=tranNp, line.length=lineLength,
                                                     edgeAdjust=FALSE, spat.random.type=rand.type, return.index=TRUE))
  
  #get lots of designs  
  combinedSamples <- parallel::mclapply( 1:nSims, FUN=sampleTransect_strat, mc.cores=mc.cores, n=n, strata=myStrata,
                                         inclusion.probs=p, potential.sites=X, study.area=studyArea, control=control)
  #combine them (matrices are apparently quicker)
  manyDesigns <- parallel::mclapply( combinedSamples, FUN=function(x) as.matrix( x))
  pointSamples <- matrix( NA, ncol=2, nrow=nSims * n * tranNp)
  for( ii in 1:nSims)
    pointSamples[(ii-1)*n*tranNp + 1:(n*tranNp),] <- manyDesigns[[ii]]
  pointSamples <- as.data.frame( pointSamples)
  
  #find the empirical (frequentist) probability of observing each cell
  pointSamples$cellID <- class::knn1( train=X, test=pointSamples, cl=1:N)
  cellVisits <- table( pointSamples$cellID) #this is already a factor with N levels.
  cellProps$strata <- cellVisits / sum( cellVisits)
  
  #calculate the MSE
  results$incl.probs.mse["mse","strata"] <- sum((cellProps$strata-p)^2)/N
  #and the abs deviation
  results$incl.probs.mse["absDev","strata"] <- mean( abs( cellProps$strata-p))/N
  
  #calculate the average spat. bal. and its sd
  voronoiVars <- parallel::mclapply( manyDesigns, getVarArea, X.grid=X, inclProbs=p, mc.cores=mc.cores)
  voronoiVars <- t( simplify2array( voronoiVars))
  
  results$spat.bal["mean_var","strata"] <- mean( voronoiVars[,"var"])
  results$spat.bal["sd_var","strata"] <- sd( voronoiVars[,"var"])
  
  tmp <- t( simplify2array( parallel::mclapply( 1:nSims, designTranDist, mc.cores=mc.cores)))
  
  results$spat.bal["mean_minTransDist","strata"] <- mean( tmp[,'min'])
  results$spat.bal["sd_minTransDist","strata"] <- sd( tmp[,'min'])
  results$spat.bal["mean_meanTransDist","strata"] <- mean( tmp[,'mean'])
  results$spat.bal["sd_meanTransDist","strata"] <- sd( tmp[,'mean'])
}


##########################################################
####	Saving results and designs
##########################################################

save( results, cellProps, file="simItems1.RData")

#########################################################
####	Plotting/Reporting Results                       
#########################################################
{
  pdf( "VolcanoSim_quasi_results2.pdf", height=7.7, width=4.45)
  
  layout( matrix( c(1,1, 4,9, 2,2, 5,7, 3,3, 6,8, 10,10), nrow=7, ncol=2, byrow=TRUE), widths=c(1.6,1), heights=c(2,5,1,5,1,5,0.5))
  
  par( mar=c( 0,0,0,0))
  plot.new()
  text(x=0.5,y=0.5,adj=c(0.5,0.25),labels="Inclusion Probabilities (BAS)", cex=2.5)
  
  plot.new()
  #  text(x=0.5,y=0.5,adj=c(0.5,1.5),labels="Observed with Naive BAS", cex=1.5)
  
  plot.new()
  #  text(x=0.5,y=0.5,adj=c(0.5,1.5),labels="Observed with Transect BAS", cex=1.5)
  
  par( mar=c( 1,1.5,1,1)+0.1)
  
  my.zlims <- range( sapply( cellProps[c(1,4)], range))
  
  image( x=10*(1:n.x)-5, y=10*(1:n.y)-5, z=matrix( p, ncol=n.y, nrow=n.x), asp=1, main="", yaxt='n', xaxt='n', xlab='', ylab='', col=tim.colors(), zlim=my.zlims)
  mtext("Specified", side=3, line=0.5, outer=FALSE, cex=1)
  mtext( "(A)", side=1, adj=-0.1, line=-1, cex=0.8)
  image( x=10*(1:n.x)-5, y=10*(1:n.y)-5, z=matrix( cellProps$cell_quasi, ncol=n.y, nrow=n.x), asp=1, main="", yaxt='n', xaxt='n', xlab='', ylab='', col=tim.colors(), zlim=my.zlims)
  mtext("Observed with Naive BAS", side=3, line=0.5, outer=FALSE, cex=1)
  mtext( "(B)", side=1, adj=-0.1, line=-1, cex=0.8)
  image( x=10*(1:n.x)-5, y=10*(1:n.y)-5, z=matrix( cellProps$Transect_quasi, ncol=n.y, nrow=n.x), asp=1, main="", yaxt='n', xaxt='n', xlab='', ylab='', col=tim.colors(), zlim=my.zlims)
  mtext("Observed with Transect BAS", side=3, line=0.5, outer=FALSE, cex=1)
  mtext( "(C)", side=1, adj=-0.1, line=-1, cex=0.8)
  
  par( mar=c( 2, 3, 1, 1)+0.1, xpd=FALSE, pty='s')
  
  combined.probs <- c(p, cellProps$cell_quasi, cellProps$Transect_quasi_noedge, cellProps$Transect_quasi)
  my.lims <- range( combined.probs)
  
  plot( p, cellProps$cell_quasi, pch=20, main="", xlab="", ylab="", asp=1, xlim=my.lims, ylim=my.lims, axes=FALSE)
  abline( 0,1, col='red')
  box()
  axis( 1, at=pretty( combined.probs, 5), labels=pretty( combined.probs, 5))
  axis( 2, at=pretty( combined.probs, 5), labels=pretty( combined.probs, 5))
  mtext("Observed", side=2, line=2, outer=FALSE, cex=0.75)
  mtext( "Specified", side=1, line=2, outer=FALSE, cex=0.75)
  mtext( "(D)", side=1, adj=-0.4, line=1.75, cex=0.8)
  
  plot( p, cellProps$Transect_quasi, pch=20, main="", xlab="", ylab="", asp=1, xlim=my.lims, ylim=my.lims, axes=FALSE)
  abline( 0,1, col='red')
  box()
  axis( 1, at=pretty( combined.probs, 5), labels=pretty( combined.probs, 5))
  axis( 2, at=pretty( combined.probs, 5), labels=pretty( combined.probs, 5))
  mtext( "Observed", side=2, line=2, outer=FALSE, cex=0.75)
  mtext( "Specified", side=1, line=2, outer=FALSE, cex=0.75)
  mtext( "(E)", side=1, adj=-0.4, line=1.75, cex=0.8)
  
  dev.off()
}

{
  pdf( "VolcanoSim_pseudo_results2.pdf", height=7.7, width=4.45)
  
  layout( matrix( c(1,1, 4,9, 2,2, 5,7, 3,3, 6,8, 10,10), nrow=7, ncol=2, byrow=TRUE), widths=c(1.6,1), heights=c(2,5,1,5,1,5,0.5))
  
  par( mar=c( 0,0,0,0))
  plot.new()
  text(x=0.5,y=0.5,adj=c(0.5,0.25),labels="Inclusion Probabilities (SRS)", cex=2.5)
  #  text(x=0.5,y=0.5,adj=c(0.5,0.25),labels="Inclusion Probabilities (SRS)", cex=2.5)
  
  plot.new()
  #  text(x=0.5,y=0.5,adj=c(0.5,1.5),labels="Observed with Naive SRS", cex=1.5)
  
  plot.new()
  #  text(x=0.5,y=0.5,adj=c(0.5,1.5),labels="Observed with Transect SRS", cex=1.5)
  
  par( mar=c( 1,1.5,1,1)+0.1)
  
  my.zlims <- range( sapply( cellProps[c(2,3)], range))
  
  image( x=10*(1:n.x)-5, y=10*(1:n.y)-5, z=matrix( p, ncol=n.y, nrow=n.x), asp=1, main="", yaxt='n', xaxt='n', xlab='', ylab='', col=tim.colors(), zlim=my.zlims)
  mtext("Specified", side=3, line=0.5, outer=FALSE, cex=1)
  mtext( "(A)", side=1, adj=-0.1, line=-1, cex=0.8)
  image( x=10*(1:n.x)-5, y=10*(1:n.y)-5, z=matrix( cellProps$cell_pseudo, ncol=n.y, nrow=n.x), asp=1, main="", yaxt='n', xaxt='n', xlab='', ylab='', col=tim.colors(), zlim=my.zlims)
  mtext("Observed with Naive SRS", side=3, line=0.5, outer=FALSE, cex=1)
  mtext( "(B)", side=1, adj=-0.1, line=-1, cex=0.8)
  image( x=10*(1:n.x)-5, y=10*(1:n.y)-5, z=matrix( cellProps$Transect_pseudo, ncol=n.y, nrow=n.x), asp=1, main="", yaxt='n', xaxt='n', xlab='', ylab='', col=tim.colors(), zlim=my.zlims)
  mtext("Observed with Transect SRS", side=3, line=0.5, outer=FALSE, cex=1)
  mtext( "(C)", side=1, adj=-0.1, line=-1, cex=0.8)
  
  par( mar=c( 2, 3, 1, 1)+0.1, xpd=FALSE, pty='s')
  
  combined.probs <- c(p, cellProps$cell_pseudo, cellProps$Transect_pseudo_noedge, cellProps$Transect_pseudo)
  my.lims <- range( combined.probs)
  
  plot( p, cellProps$cell_pseudo, pch=20, main="", xlab="", ylab="", asp=1, xlim=my.lims, ylim=my.lims, axes=FALSE)
  abline( 0,1, col='red')
  box()
  axis( 1, at=pretty( combined.probs, 5), labels=pretty( combined.probs, 5))
  axis( 2, at=pretty( combined.probs, 5), labels=pretty( combined.probs, 5))
  mtext("Observed", side=2, line=2, outer=FALSE, cex=0.75)
  mtext( "Specified", side=1, line=2, outer=FALSE, cex=0.75)
  mtext( "(D)", side=1, adj=-0.4, line=1.75, cex=0.8)
  
  plot( p, cellProps$Transect_pseudo, pch=20, main="", xlab="", ylab="", asp=1, xlim=my.lims, ylim=my.lims, axes=FALSE)
  abline( 0,1, col='red')
  box()
  axis( 1, at=pretty( combined.probs, 5), labels=pretty( combined.probs, 5))
  axis( 2, at=pretty( combined.probs, 5), labels=pretty( combined.probs, 5))
  mtext( "Observed", side=2, line=2, outer=FALSE, cex=0.75)
  mtext( "Specified", side=1, line=2, outer=FALSE, cex=0.75)
  mtext( "(E)", side=1, adj=-0.4, line=1.75, cex=0.8)
  
  dev.off()
}

{
  pdf( "VolcanoSim_strata_results2.pdf", height=7.7, width=4.45)
  
  layout( matrix( c(1,1, 4,9, 2,2, 5,7, 3,3, 6,8, 10,10), nrow=7, ncol=2, byrow=TRUE), widths=c(1.6,1), heights=c(2,5,1,5,1,5,0.5))
  
  par( mar=c( 0,0,0,0))
  plot.new()
  text(x=0.5,y=0.5,adj=c(0.5,0.25),labels="Inclusion Probabilities (stratified)", cex=2.5)
  #  text(x=0.5,y=0.5,adj=c(0.5,0.25),labels="Inclusion Probabilities (SRS)", cex=2.5)
  
  plot.new()
  #  text(x=0.5,y=0.5,adj=c(0.5,1.5),labels="Observed with Naive SRS", cex=1.5)
  
  plot.new()
  #  text(x=0.5,y=0.5,adj=c(0.5,1.5),labels="Observed with Transect SRS", cex=1.5)
  
  par( mar=c( 1,1.5,1,1)+0.1)
  
  my.zlims <- range( sapply( cellProps[c(2,3)], range))
  
  image( x=10*(1:n.x)-5, y=10*(1:n.y)-5, z=matrix( p, ncol=n.y, nrow=n.x), asp=1, main="", yaxt='n', xaxt='n', xlab='', ylab='', col=tim.colors(), zlim=my.zlims)
  mtext("Specified", side=3, line=0.5, outer=FALSE, cex=1)
  mtext( "(A)", side=1, adj=-0.1, line=-1, cex=0.8)
  strata <- findStrata_fromInclProbs( x=p, nstrata = nstrata)
  image( x=10*(1:n.x)-5, y=10*(1:n.y)-5, z=matrix( strata, ncol=n.y, nrow=n.x), asp=1, main="", yaxt='n', xaxt='n', xlab='', ylab='', col=tim.colors())
  mtext("Strata", side=3, line=0.5, outer=FALSE, cex=1)
  mtext( "(B)", side=1, adj=-0.1, line=-1, cex=0.8)
  image( x=10*(1:n.x)-5, y=10*(1:n.y)-5, z=matrix( cellProps$strata, ncol=n.y, nrow=n.x), asp=1, main="", yaxt='n', xaxt='n', xlab='', ylab='', col=tim.colors(), zlim=my.zlims)
  mtext("Observed with Stratified", side=3, line=0.5, outer=FALSE, cex=1)
  mtext( "(C)", side=1, adj=-0.1, line=-1, cex=0.8)
  
  plot.new()
  
  par( mar=c( 2, 3, 1, 1)+0.1, xpd=FALSE, pty='s')
  
  combined.probs <- c(p, cellProps$strata)
  my.lims <- range( combined.probs)
  
  plot( p, cellProps$strata, pch=20, main="", xlab="", ylab="", asp=1, xlim=my.lims, ylim=my.lims, axes=FALSE)
  abline( 0,1, col='red')
  box()
  axis( 1, at=pretty( combined.probs, 5), labels=pretty( combined.probs, 5))
  axis( 2, at=pretty( combined.probs, 5), labels=pretty( combined.probs, 5))
  mtext("Observed", side=2, line=2, outer=FALSE, cex=0.75)
  mtext( "Specified", side=1, line=2, outer=FALSE, cex=0.75)
  mtext( "(D)", side=1, adj=-0.4, line=1.75, cex=0.8)
  
  dev.off()
}