
############### Examples of MBH design package from CRAN PDF ####################

rm( list=ls())
#.libPaths( c(.libPaths(), "~/lib/R/library"))

library( MBHdesign)
library( parallel)
library( class)
library( fields)
library( pdist)


set.seed( 747) #I'm currently on a 787, so it *almost* seems appropriate
#number of transects
n <- 10
#number of points to sample from
N <- 100^2
#the sampling grid (offset so that the edge locations have same area)
offsetX <- 1/(2*sqrt( N))
my.seq <- seq( from=offsetX, to=1-offsetX, length=sqrt(N))
X <- expand.grid( my.seq, my.seq)
colnames( X) <- c("X1","X2")

#############
## Inclusion probabilities -----
# Pretend there is some gradi#non-uniform inclusion probabilities
inclProbs <- 1-exp(-X[,1])
#scaling to enforce summation to n
inclProbs <- n * inclProbs / sum( inclProbs)
#uniform inclusion probabilities would be inclProbs <- rep( n/N, times=N)
#visualise
image( x=unique( X[,1]), y=unique( X[,2]),
       z=matrix( inclProbs, nrow=sqrt(nrow(X)), ncol=sqrt(nrow( X))),
       main="(Undadjusted) Inclusion Probabilities",
       ylab=colnames( X)[2], xlab=colnames( X)[1])

####################
## Transect sampling funciton ----
# Transects as a set of points in a line 
# This example is in favour of computation, not accuracy
# there are other control options about the algorithm

#my.control is a list that contains
my.control <- list(
  #the type of transect
  transect.pattern="line",
  #the length of transect
  line.length=0.15,
  #the number of points that define the transect
  transect.nPts=15,
  #the number of putative directions that a transect can take
  nRotate=9
)

# take the transect sample --
samp <- transectSamp( n=n, potential.sites=X, inclusion.probs=inclProbs,
                      control=my.control)
image( x=unique( X[,1]), y=unique( X[,2]),
       z=matrix( inclProbs, nrow=sqrt(nrow(X)), ncol=sqrt(nrow( X))),
       main="(Undadjusted) Inclusion Probabilities",
       sub="10 Transects",
       ylab=colnames( X)[2], xlab=colnames( X)[1])
points( samp$points[,5:6], pch=20, cex=0.6)

## These can now be written as a csv file --
#write csv
#write.csv( samp$transect, file="transectSample1.csv", row.names=FALSE)
#tidy
#rm( list=ls())


