# ------------------------------------------------
# DensityAroundOnePoint.R
# ------------------------------------------------
# Code provided as is and can be used or modified freely. 
#
# We would however appreciate if you cite the following paper:
#
# Demsar U, Buchin K, van Loon EE and Shamoun-Baranes J, 2015, 
# Stacked space-time densities: a geovisualisation approach to explore 
# dynamics of space use over time. GeoInformatica, 19(1):85-115.
# DOI 10.1007/s10707-014-0207-5
# ------------------------------------------------
# Author: Urska Demsar
# University of St Andrews
# St Andrews, Scotland, UK
# http://udemsar.com
# Date: 21 March 2016
# ------------------------------------------------
# Supporting file to StackedSpaceTimeDensities.R
# ------------------------------------------------
# Parameters:
# currentLine = one point as a row of 3D coordinates 
# kernelSize = size of kernel for distance cut off, in original units
# xcoord, ycoord, zcoord = 3D raster of coordinates in each direction
# method = type of kernel: 1 = linear, 2 = bisquare, 3 = Gaussian, 4 = Brownian bridge
# voxelSize = size of voxel (resolution)
# sigma11 = variance 1 for Brownian bridges (movement uncertainty)
# sigma12 = variance 2 for Brownian bridges (sensor location uncertainty)
# ------------------------------------------------
# Output:
# vcoord - volume with density values around the one point - only in the layer where the point is
# ------------------------------------------------

DensityAroundOnePoint <- function(currentLine,kernelSize,xcoord,ycoord,zcoord,method,voxelSize,sigma11,sigma12) {

  #head(currentLine)
  
  # read the size of the 3D rasters - xcoord, ycoord and zcoord must be of 
  # same size, plus calculate the range in each direction.
  sizeOf3Draster <- dim(xcoord)
  
  xnn <- sizeOf3Draster[1]-1
  ynn <- sizeOf3Draster[2]-1
  znn <- sizeOf3Draster[3]-1
  
  minXcoord <- min(xcoord[,1,1])
  maxXcoord <- max(xcoord[,1,1])
  Xrange <- maxXcoord-minXcoord
  
  minYcoord <- min(ycoord[1,,1])
  maxYcoord <- max(ycoord[1,,1])
  Yrange <- maxYcoord-minYcoord
  
  minZcoord <- min(zcoord[1,1,])
  maxZcoord <- max(zcoord[1,1,])
  Zrange <- maxZcoord-minZcoord
  
  kx <- xnn/Xrange
  ky <- ynn/Yrange
  kz <- znn/Zrange
  
  nx <- (maxXcoord-(xnn+1)*minXcoord)/Xrange
  ny <- (maxYcoord-(ynn+1)*minYcoord)/Yrange
  nz <- (maxZcoord-(znn+1)*minZcoord)/Zrange
  
  # initialise an empty vcoord (this will eventually be density volume):
  vcoord <- array(data=0,dim=dim(xcoord))
  
  # Calculate the density by finding the layer of the given point and calculate density around it
  
  # Take the point P(xn,yn,zn) 
  xn <- currentLine[1,1]
  yn <- currentLine[1,2]
  zn <- currentLine[1,3]
    
  # Find the temporal voxel layer of this point
  
  k <- round(kz*zn+nz)
  
  # Calculate parameters for BBs - adapted from the segment version, to match the parameters when the trajectory
  # has more than two points
  if (method==4) {
    # for the segment it was: TT <- zn1-zn, this can be at min the size of a voxel
    TT <- voxelSize
    p <- 0.5 # this was interpolation between P and P1, now we set it to the centre of vozel layer
    sigma2 <- sigma11*(1-p)*p*TT + sigma12*TT*((1-p)^2+p^2)
  }
  
  # Calculate density in this one layer around the point
  for (i in 1:xnn+1) {
    for (j in 1:ynn+1) {
          
      # Calculate the distance from each voxel in the bounding box 
      # around the line segment to the interpolated point on the
      # line segment at that temporal level (i.e. at that k)
      # (x,y,z)= voxel point (i,j,k)
      # (xn,yn,zn) = original point at level k
        x <- xcoord[i,j,k]
        y <- ycoord[i,j,k]
        z <- zcoord[i,j,k]
        dist <- sqrt(SquaredDistance2points(x,y,z,xn,yn,zn))
      
      # Linear decay function
          if (method==1) {
            if (dist>kernelSize) {
              vcoord[i,j,k]=0 
            }
            else {
              vcoord[i,j,k]=LinearDecay(dist,kernelSize,voxelSize) 
            }
          } # if (method==1)
          
      # Bisquare decay function
          if (method==2) {
            if (dist>kernelSize) {
              vcoord[i,j,k]=0
            }
            else {
              vcoord[i,j,k]=BisquareDecay(dist,kernelSize,voxelSize)
            }
          } # if (method==2)
          
      # Gaussian decay function
          # For visual comparison with other densities I choose
          # mu and sigma2 to be such that the value at distance
          # 0 is 1, and at kernel distance very small (I chose
          # exp(-1). This is at mu = 0 and sigma2=(kernelSize/2)^2
          if (method==3) { 
            vcoord[i,j,k]=GaussianDecay(dist,0,(kernelSize/2)^2,voxelSize)
          } # if (method==3)
          
      # Brownian bridges decay function
          # For Brownian decay we have a different Gaussian
          # density at each level, the width of the Gaussian
          # density depends on sigma2 (see Bullard 1991 and Horne
          # et al. 2007 for details on calculation)
          # If there is only one point, we take the original largest Gaussian 
          if (method==4) {
            vcoord[i,j,k]=GaussianDecay(dist,0,sigma2,voxelSize)     
          } # if (method==4)
           
        } # for (j in 1:ynn+1)
      } # for (i in 1:xnn+1)
      
# Testing the level of non-zero layer, if it's correct, it should be the same as k
# nonZeroVcoord <- which(vcoord > 0)
# nonZeroZcoord <- zcoord[nonZeroVcoord]/100
# unique(nonZeroZcoord)-k
# OK! Let's also check it in voxler -> export
#  WriteCSVtableFromFourRectangular3DArrays(xcoord,ycoord,zcoord,vcoord,'ID_4_OneDayOnly_1160_OnePointInADay_Density.csv')
  
return(vcoord)
}