#############################################################################################
####	Visual illustration for the steps in designing transects
####	Final results were produced using R-3.6.1 during July 2019
####	We note that the random number generator was updated for 
####	R-3.6.x and so results may differ if script is run on a previous
####	version of R.
####################################################################

jet.seed <- 747 #seed for randomisation

library( MBHdesign)
library( fields)

#generate samples on a grid
#the number of potential sampling locations
N <- 100^2
#number of samples
n <- 10
#the grid on unit square
X <- as.matrix( expand.grid( 1:sqrt( N), 1:sqrt(N)) / sqrt(N) - 1/(2*sqrt(N)))
##the inclusion probabiltiies with gradient according to non-linear function of X[,1]
p <- dcauchy( X[,1], location=0.66, scale=0.5) * dcauchy( X[,2], location=0.66, scale=0.5)
p <- pnorm( ( X[,1] + X[,2])) * p
#standardise to get n samples
p <- p / sum( p)

#get the sample

nRotate <- 51
transect.nPts <- 20
transLength <- 0.25

set.seed( jet.seed)

#code ripped from transectSamp for this purpose
control <- MBHdesign:::set.transect.control( list( nRotate=nRotate, transect.nPts=transect.nPts, line.length=transLength, 
                                                   edgeAdjust=TRUE, mc.cores=8))
#not calling id.list -- nothing to get. Specify it instead
id.list <- list( IDs=NULL, adjustedSpecified=NULL)

#make the cluster for computing
cl <- parallel::makeCluster(control$mc.cores)
parallel::clusterExport(cl, c( "getSingleTransLocProbs", "getSingleSiteLocProbs"), envir=.getNamespace("MBHdesign"))
parallel::clusterExport( cl, "in.out", envir=.getNamespace("mgcv"))
cl.flag <- TRUE

#set the parameters of the design (survey area, potential sites, inclusion probs, all the different n's, ...)
designParams <- MBHdesign:::setDesignParams_transect( study.area=matrix( c( 0,0, 1,0, 1,1, 0,1), ncol=2, byrow=TRUE), potential.sites=X, inclusion.probs=p, control=control)
#set parameters of the individual transects (pattern, direction, point representation, ...)
transectParams <- MBHdesign:::setTransectParams2( designParams=designParams, control=control)
#cell/transect inclusion probs and transect composition (which cells in which transect)
locProbs.raw <- MBHdesign:::getProbs2( designParams, transectParams, constraints=NULL, mc.cores=control$mc.cores, prev.ids=id.list$IDs, cl=cl)

#Edge adjusted inclusion probabilities for transects (nSites X nRotates) and cells (nSites X 1)
locProbs.edge <- MBHdesign:::adjustEdge4( locProbs.raw, designParams, transectParams, constrainedSet=NULL, control=control, cl=cl, 
                                          adjustedSpecified=id.list$adjustedSpecified)
control$edgeAdjust <- FALSE
control$calcObsProbs <- TRUE
locProbs.edge.noAdjust <- MBHdesign:::adjustEdge4( locProbs.raw, designParams, transectParams, constrainedSet=NULL, control=control)

#plotting probabilities

pdf( "~/NESP/MonitoringTheme/SOPproject/SpatDesignStuff/TransectDesign/Simulations/VisualIllustration/visualisation5.pdf", height=7, width=8.75)

layout( matrix( c( 1,1,1, 2,3,4, 5,6,7), ncol=3, nrow=3, byrow=TRUE), widths=c(1,1,1,0.4), heights=c(1,5,5))

par( mar=c( 0,0,0,0))
plot.new()
text(x=0.5,y=0.5,adj=c(0.5,0.25),labels="Transect Sampling Process", cex=2.5)

#my.zlims <- range( c(p, rowSums( locProbs.raw$transect.probs), rowSums( locProbs.edge$transects)))
par( mar=c( 2,2,2,0)+0.1)
image( x=designParams$potential.sites.list$Var1, y=designParams$potential.sites.list$Var2, 
       z=matrix( p, nrow=designParams$N$dimension1, ncol=designParams$N$dimension2),
       xlab="", ylab="", asp=1, bty='o', col=tim.colors(), pty='s', xaxt='n', yaxt='n', main="Specified")
mtext( "(A)", side=1, adj=-0.1, line=-1, cex=0.8)

image( x=designParams$potential.sites.list$Var1, y=designParams$potential.sites.list$Var2, 
       z=matrix( attr( locProbs.edge.noAdjust, "observed"), nrow=designParams$N$dimension1, ncol=designParams$N$dimension2),
       main="", xlab="", ylab="", asp=1, bty='o', col=tim.colors(), pty='s', xaxt='n', yaxt='n')
title(main="Observed (No Adjustment)", line=1)
mtext( "(B)", side=1, adj=-0.1, line=-1, cex=0.8)

image( x=designParams$potential.sites.list$Var1, y=designParams$potential.sites.list$Var2, 
       z=matrix( attr( locProbs.edge, "observed"), nrow=designParams$N$dimension1, ncol=designParams$N$dimension2),
       main="", xlab="", ylab="", asp=1, yaxt='n', bty='o', col=tim.colors(), pty='s', xaxt='n')
title(main="Observed (Edge Adjusted)", line=1)
mtext( "(C)", side=1, adj=-0.1, line=-1, cex=0.8)

image( x=designParams$potential.sites.list$Var1, y=designParams$potential.sites.list$Var2, 
       z=matrix( rowSums( locProbs.edge$transects, na.rm=TRUE), nrow=designParams$N$dimension1, ncol=designParams$N$dimension2),
       xlab="", ylab="", asp=1, bty='o', col=tim.colors(), pty='s', xaxt='n', yaxt='n', main="Edge-Adjusted Transect")
mtext( "(D)", side=1, adj=-0.1, line=-1, cex=0.8)

#take the sample and massage it into form
transIDs <- MBHdesign:::getTransects_2step( n, locProbs.edge$transect, designParams, transectParams, control)
ret <- list()
ret$transect <- cbind( 1:n, do.call( "rbind", lapply( transIDs, function(x) c(x$start_locat, x$direction, x$inclusion.prob.transect*n))))
colnames( ret$transect) <- c("transect", names( transIDs[[1]]$start_locat), "direction", "transect_inclusionProb")
ret$points <- lapply( transIDs, function(x) x$locats)
ret$points <- lapply( 1:length( ret$points), function(x) cbind( x, ret$transect[x, names( transIDs[[1]]$start_locat)[1]], ret$transect[x, names( transIDs[[1]]$start_locat)[2]], ret$transect[x, "direction"], ret$points[[x]]))
ret$points <- do.call( "rbind", ret$points)
ret$points <- cbind( ret$points, do.call( 'c', lapply( transIDs, function(x) locProbs.raw$cell.probs[x$site, x$rotate, ])))
ret$points <- cbind( ret$points, locProbs.edge$cells[do.call( 'c', lapply( transIDs, function(x) locProbs.raw$IDs[x$site, x$rotate, ]))])
colnames( ret$points) <- c( "transect", paste0( "start_", names( transIDs[[1]]$start_locat)), "direction", names( transIDs[[1]]$start_locat), "specifiedInclProb", "AdjustedInclProb")
ret$points <- as.data.frame( ret$points)
ret$transect <- as.data.frame( ret$transect)

#par( mar=c( 2,2,2,0)+0.1)
image( x=(1:sqrt(N))/sqrt(N), y=(1:sqrt(N))/sqrt(N), z=matrix( p, ncol=sqrt(N), nrow=sqrt(N)), 
       main="Transect Centres", xlab="", ylab="", asp=1, bty='o', col=tim.colors(), pty='s', xaxt='n', yaxt='n')
points( ret$transect[,c("Var1","Var2")], cex=1.5)
mtext( "(E)", side=1, adj=-0.1, line=-1, cex=0.8)

image( x=(1:sqrt(N))/sqrt(N), y=(1:sqrt(N))/sqrt(N), z=matrix( p, ncol=sqrt(N), nrow=sqrt(N)), 
       main="All Transects", xlab="", ylab="", asp=1, bty='o', col=tim.colors(), pty='s', xaxt='n', yaxt='n')
points( ret$points[,c("Var1","Var2")], pch=20, cex=1.5)
points( ret$transect[,c("Var1","Var2")], cex=1.5)
mtext( "(F)", side=1, adj=-0.1, line=-1, cex=0.8)

dev.off()


###########################################################################################
####  An 's' transect shape on unequal inclusion prob gird.

pattern <- matrix( c(0,0, 0.05,0, 0.05,0.05, 0,0.05, 0,0.1, 0.05,0.1), ncol=2, byrow=TRUE)
m <- 4
pattern1 <- matrix( NA, ncol=2, nrow=(nrow( pattern)-1)*(m))
for( ii in 1:(nrow( pattern)-1)){
  pattern1[(ii-1)*m+1:m,1] <- seq( from=pattern[ii,1], to=pattern[ii+1,1], length=m)
  pattern1[(ii-1)*m+1:m,2] <- seq( from=pattern[ii,2], to=pattern[ii+1,2], length=m)
}
pattern1 <- pattern1[!duplicated( pattern1),]

set.seed(jet.seed)

#code ripped from transectSamp for this purpose
control <- MBHdesign:::set.transect.control( list( nRotate=nRotate, transect.nPts=transect.nPts, line.length=transLength, 
                                                   edgeAdjust=TRUE, mc.cores=8))
#not calling id.list -- nothing to get. Specify it instead
id.list <- list( IDs=NULL, adjustedSpecified=NULL)

#code ripped from transectSamp for this purpose
control <- MBHdesign:::set.transect.control( list( nRotate=nRotate, transect.nPts=transect.nPts, line.length=0.1, edgeAdjust=TRUE, mc.cores=8, transect.pattern=pattern1))
#not calling id.list -- nothing to get. Specify it instead
id.list <- list( IDs=NULL, adjustedSpecified=NULL)

#make the cluster for computing
cl <- parallel::makeCluster(control$mc.cores)
parallel::clusterExport(cl, c( "getSingleTransLocProbs", "getSingleSiteLocProbs"), envir=.getNamespace("MBHdesign"))
parallel::clusterExport( cl, "in.out", envir=.getNamespace("mgcv"))
cl.flag <- TRUE

#set the parameters of the design (survey area, potential sites, inclusion probs, all the different n's, ...)
designParams <- MBHdesign:::setDesignParams_transect( study.area=matrix( c( 0,0, 1,0, 1,1, 0,1), ncol=2, byrow=TRUE), potential.sites=X, inclusion.probs=p, control=control)
#set parameters of the individual transects (pattern, direction, point representation, ...)
transectParams <- MBHdesign:::setTransectParams2( designParams=designParams, control=control)
#cell/transect inclusion probs and transect composition (which cells in which transect)
locProbs.raw <- MBHdesign:::getProbs2( designParams, transectParams, constraints=NULL, mc.cores=control$mc.cores, prev.ids=id.list$IDs, cl=cl)

#Edge adjusted inclusion probabilities for transects (nSites X nRotates) and cells (nSites X 1)
locProbs.edge <- MBHdesign:::adjustEdge4( locProbs.raw, designParams, transectParams, constrainedSet=NULL, control=control, cl=cl, adjustedSpecified=id.list$adjustedSpecified)
control$edgeAdjust <- FALSE
control$calcObsProbs <- TRUE
locProbs.edge.noAdjust <- MBHdesign:::adjustEdge4( locProbs.raw, designParams, transectParams, constrainedSet=NULL, control=control)

#plotting probabilities

pdf( "~/NESP/MonitoringTheme/SOPproject/SpatDesignStuff/TransectDesign/Simulations/VisualIllustration/visualisation_shape5.pdf", height=7, width=8.75)

layout( matrix( c( 1,1,1, 2,3,4, 5,6,7), ncol=3, nrow=3, byrow=TRUE), widths=c(1,1,1,0.4), heights=c(1,5,5))

par( mar=c( 0,0,0,0))
plot.new()
text(x=0.5,y=0.5,adj=c(0.5,0.25),labels="Transect Sampling Process", cex=2.5)

#my.zlims <- range( c(p, rowSums( locProbs.raw$transect.probs), rowSums( locProbs.edge$transects)))
par( mar=c( 2,2,2,0)+0.1)
image( x=designParams$potential.sites.list$Var1, y=designParams$potential.sites.list$Var2, 
       z=matrix( p, nrow=designParams$N$dimension1, ncol=designParams$N$dimension2),
       xlab="", ylab="", asp=1, bty='o', col=tim.colors(), pty='s', xaxt='n', yaxt='n', main="Specified")
mtext( "(A)", side=1, adj=-0.1, line=-1, cex=0.8)

image( x=designParams$potential.sites.list$Var1, y=designParams$potential.sites.list$Var2, 
       z=matrix( attr( locProbs.edge.noAdjust, "observed"), nrow=designParams$N$dimension1, ncol=designParams$N$dimension2),
       main="", xlab="", ylab="", asp=1, bty='o', col=tim.colors(), pty='s', xaxt='n', yaxt='n')
title(main="Observed (No Adjustment)", line=1)
mtext( "(B)", side=1, adj=-0.1, line=-1, cex=0.8)

image( x=designParams$potential.sites.list$Var1, y=designParams$potential.sites.list$Var2, 
       z=matrix( attr( locProbs.edge, "observed"), nrow=designParams$N$dimension1, ncol=designParams$N$dimension2),
       main="", xlab="", ylab="", asp=1, yaxt='n', bty='o', col=tim.colors(), pty='s', xaxt='n')
title(main="Observed (Edge Adjusted)", line=1)
mtext( "(C)", side=1, adj=-0.1, line=-1, cex=0.8)

image( x=designParams$potential.sites.list$Var1, y=designParams$potential.sites.list$Var2, 
       z=matrix( rowSums( locProbs.edge$transects, na.rm=TRUE), nrow=designParams$N$dimension1, ncol=designParams$N$dimension2),
       xlab="", ylab="", asp=1, bty='o', col=tim.colors(), pty='s', xaxt='n', yaxt='n', main="Edge-Adjusted Transect")
mtext( "(D)", side=1, adj=-0.1, line=-1, cex=0.8)

#take the sample and massage it into form
transIDs <- MBHdesign:::getTransects_2step( n, locProbs.edge$transect, designParams, transectParams, control)
ret <- list()
ret$transect <- cbind( 1:n, do.call( "rbind", lapply( transIDs, function(x) c(x$start_locat, x$direction, x$inclusion.prob.transect*n))))
colnames( ret$transect) <- c("transect", names( transIDs[[1]]$start_locat), "direction", "transect_inclusionProb")
ret$points <- lapply( transIDs, function(x) x$locats)
ret$points <- lapply( 1:length( ret$points), function(x) cbind( x, ret$transect[x, names( transIDs[[1]]$start_locat)[1]], ret$transect[x, names( transIDs[[1]]$start_locat)[2]], ret$transect[x, "direction"], ret$points[[x]]))
ret$points <- do.call( "rbind", ret$points)
ret$points <- cbind( ret$points, do.call( 'c', lapply( transIDs, function(x) locProbs.raw$cell.probs[x$site, x$rotate, ])))
ret$points <- cbind( ret$points, locProbs.edge$cells[do.call( 'c', lapply( transIDs, function(x) locProbs.raw$IDs[x$site, x$rotate, ]))])
colnames( ret$points) <- c( "transect", paste0( "start_", names( transIDs[[1]]$start_locat)), "direction", names( transIDs[[1]]$start_locat), "specifiedInclProb", "AdjustedInclProb")
ret$points <- as.data.frame( ret$points)
ret$transect <- as.data.frame( ret$transect)

#par( mar=c( 2,2,2,0)+0.1)
image( x=(1:sqrt(N))/sqrt(N), y=(1:sqrt(N))/sqrt(N), z=matrix( p, ncol=sqrt(N), nrow=sqrt(N)), 
       main="Transect Centres", xlab="", ylab="", asp=1, bty='o', col=tim.colors(), pty='s', xaxt='n', yaxt='n')
points( ret$transect[,c("Var1","Var2")], cex=1.5, pch=20)
mtext( "(E)", side=1, adj=-0.1, line=-1, cex=0.8)

image( x=(1:sqrt(N))/sqrt(N), y=(1:sqrt(N))/sqrt(N), z=matrix( p, ncol=sqrt(N), nrow=sqrt(N)), 
       main="All Transects", xlab="", ylab="", asp=1, bty='o', col=tim.colors(), pty='s', xaxt='n', yaxt='n')
points( ret$points[,c("Var1","Var2")], pch=20, cex=1.5)
mtext( "(F)", side=1, adj=-0.1, line=-1, cex=0.8)

dev.off()

rm( list=ls())

