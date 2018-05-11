# ------------------------------------------------
# DensityAroundTrajectory.R
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
# Date: 9 March 2015
# ------------------------------------------------
# Supporting file to StackedSpaceTimeDensities.R
# ------------------------------------------------
# Parameters:
# currentLine = trajectory as a matrix of 3D points (each row = 1 point)
# kernelSize = size of kernel for distance cut off, in original units
# xcoord, ycoord, zcoord = 3D raster of coordinates in each direction
# method = type of kernel: 1 = linear, 2 = bisquare, 3 = Gaussian, 4 = Brownian bridge
# voxelSize = size of voxel (resolution)
# sigma11 = variance 1 for Brownian bridges (movement uncertainty)
# sigma12 = variance 2 for Brownian bridges (sensor location uncertainty)
# ------------------------------------------------
# Output:
# vcoord - volume with density values around the trajectory
# ------------------------------------------------

DensityAroundTrajectory <- function(currentLine,kernelSize,xcoord,ycoord,zcoord,method,voxelSize,sigma11,sigma12) {

  head(currentLine)
  
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
  
  # get the length of the trajectory
  lll <- dim(currentLine)
  points <- lll[1]
  col <- lll[2]
  
  # Calculate the density by going through the trajectory segment by segment
  
  n <- 1
  
  while (n < points) {
    
    # Subsample the trajectory to voxel rows: 
    # Take two consecutive points in the trajectory, P(xn,yn,zn) and P1(xn1,yn1,zn1)
    # where P(xn1,yn1,zn1) is the first point that is temporaly away for more
    # than a voxel size than point P1(xn,yn,zn), i.e. zn1-zn>deltaZ
    xn <- currentLine[n,1]
    yn <- currentLine[n,2]
    zn <- currentLine[n,3]
    
    # Find the first next trajectory point that is in the next voxel row.
    n1 <- n+1
    while ((n1<points) && (currentLine[n1,3]-currentLine[n,3])<=voxelSize) {
      n1 <- n1+1
    }
    xn1 <- currentLine[n1,1]
    yn1 <- currentLine[n1,2]
    zn1 <- currentLine[n1,3]
    
    # Define the temporal boundaries around the segment within which the density
    # will be calculated
    kn <- min(round(kz*zn+nz),round(kz*zn1+nz))
    kn1 <- max(round(kz*zn+nz),round(kz*zn1+nz))
    
    # Now loop through all voxel layers that intersect the segment P-P1
    for (k in kn:kn1) {
      
      # Interpolate the point at level k = interpolated between P and
      # P1, this point has coordinates (xc, yc, zc).
      # p = interpolation parameter on the line segment, [0,1] from P to P1
      # If however kn and kn1 are in the same voxel layer, then just take the
      # middle point between P and P1 and don't interpolate.
      if (kn!=kn1) {
      p <- k/(kn1-kn)-kn/(kn1-kn)
      }
      else {p <- 0.5}
      
      xc <- xn+p*(xn1-xn)
      yc <- yn+p*(yn1-yn)
      zc <- zn+p*(zn1-zn)
      
      # Calculate additional parameters for Brownian density:
      # TT = time between P and P1
      # sigma2 = describes contributions of the two variances, sigma12 and
      # sigma12 to the general variance at interpolated point (given via
      # parameter p, where movement went over the entire time T), where:
      # sigma11 = variance depending on the speed of the movement
      # sigma12 = variance describing the uncertainty around end points
      if (method==4) {
        TT <- zn1-zn
        sigma2 <- sigma11*(1-p)*p*TT + sigma12*TT*((1-p)^2+p^2)
      }

      # Calculate density
      for (i in 1:xnn+1) {
        for (j in 1:ynn+1) {
          
          # Calculate the distance from each voxel in the bounding box 
          # around the line segment to the interpolated point on the
          # line segment at that temporal level (i.e. at that k)
          # (x,y,z)= voxel point (i,j,k)
          # (xc,yc,zc) = interpolated point on line segment at level k
          x <- xcoord[i,j,k]
          y <- ycoord[i,j,k]
          z <- zcoord[i,j,k]
          dist <- sqrt(SquaredDistance2points(x,y,z,xc,yc,zc))
          
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
          if (method==4) {      
            vcoord[i,j,k]=GaussianDecay(dist,0,sigma2,voxelSize)     
          } # if (method==4)
           
        } # for (j in 1:ynn+1)
      } # for (i in 1:xnn+1)
      
      
    } # for (k in kn:kn1)
    
    # Move to the next segment in trajectory
    n <- n1 
    
  } # while n < points
  
return(vcoord)
}