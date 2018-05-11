# ------------------------------------------------
# Stacked space-time densities around a set of trajectories
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
# Date: 22 March 2016
# ------------------------------------------------
# This is the main file that implements the space-time density calculation
# ------------------------------------------------
# Supporting R files:
# DecayFunctions.R
# DensityAroundOnePoint.R
# DensityAroundTrajectory.R
# SquaredDistance2points.R
# WriteCSVtableFromRectangular3DArrays.R
# ------------------------------------------------
# Input: data with a list of trajectories
# This must be a csv file including the following five variables:
# ID, ID_traj, xcoord, ycoord, zcoord, where
# - ID - ID of moving object
# - ID_traj - ID of each trajectory
# - X - x coordinate of each trajectory point in some projected
# coordinate system (in metres)
# - Y - y coordinate of each trajectory point in some projected
# coordinate system (in metres)
# - T - time of the day  
# ------------------------------------------------
# The following also need to be manually specified: 
# - volume extent (min/max X, min/max Y, min/max T)
# - voxel size
# - kernel size 
# - kernel type
# ------------------------------------------------
# Output: a csv file with density given as volume, i.e. with four columns:
# x, y, t, density value
# ------------------------------------------------


# ------------------------------------------------
# Parameters to be set/changed by the user

# Set working directory
setwd('D:/Research/3Ddensity_Code/R')

# Set voxel size in metres 
voxelSize <- 0.2

# Set kernel size in metres
kernelSize <- 1

# Set spatio-temporal extent of the volume
# This is done manually and not by reading the min/max coordinates from the
# trajectories file, so the program can be used for densities calculated
# for e.g. different days, where each day the trajectories have different
# spatial extents, but we want to add them up at the end into one general 
# space-time cube in which all density volumes must have the same extent.
#
# X & Y axis - extent in metres:
minXcoord <- 0 # westernmost point
maxXcoord <- 10 # easternmost point
minYcoord <- 0 # southernmost point
maxYcoord <- 10 # northernmost point
# Z axis - temporal extent, in seconds or minutes
minZcoord <- 0 # start point in time
maxZcoord <- 10 # end point in time

# Select kernel type: 1-linear, 2-bisquare, 3-Gaussian, 4-Brownian bridges
#method <- 1
#method <- 2
#method <- 3
method <- 4

# If Brownian bridges, the user sets sigmas 1 and 2
if (method == 4) {
  #sigma11 <- 1000 # for objects
  #sigma12 <- 300 # for objects
  sigma11 <- 3 # to be set by the user - this is just for testing
  sigma12 <- 1 # to be set by the user - this is just for testing
} else {
  sigma11 <- 0
  sigma12 <- 0
}

# Name of input and output files
# Input
trajectories.file <- 'TestData_SpatioTemporalHotspot.csv' 
# Output
#output.file <- 'TestDensity_Linear.csv'
#output.file <- 'TestDensity.Bisquare.csv'
#output.file <- 'TestDensity_Gaussian.csv'
output.file <- 'TestDensity_Brownian.csv'

# ------------------------------------------------
# General setup & read data

# Read fuctions in additional files
# Function to calculate density for each trajectory 
source("DensityAroundOnePoint.R") 
source("DensityAroundTrajectory.R") 
# Functions for saving density as as csv coordinate file
source("WriteCSVtableFromRectangular3DArrays.R") 
# Functions for geometric and kernel calculations
source("SquaredDistance2points.R")
source("DecayFunctions.R")

# Package for 3D plotting & colour schemes
library(plot3D)
library(RColorBrewer)

# Read data
dfAll <- read.csv(trajectories.file,stringsAsFactors=FALSE)
head(dfAll)

# ------------------------------------------------
# Build three 3D arrays of x, y, z coordinates and initialise the total density volume

startX <- floor(minXcoord/voxelSize)*voxelSize
startY <- floor(minYcoord/voxelSize)*voxelSize
startZ <- floor(minZcoord/voxelSize)*voxelSize
endX <- ceiling(maxXcoord/voxelSize)*voxelSize
endY <- ceiling(maxYcoord/voxelSize)*voxelSize
endZ <- ceiling(maxZcoord/voxelSize)*voxelSize

xvoxels <- (endX-startX)/voxelSize
yvoxels <- (endY-startY)/voxelSize
zvoxels <- (endZ-startZ)/voxelSize

# Build the volume
x <- seq(startX,endX,by=voxelSize)
y <- seq(startY,endY,by=voxelSize)
z <- seq(startZ,endZ,by=voxelSize)
M <- mesh(x,y,z)
xcoord <- M$x
ycoord <- M$y
zcoord <- M$z

# Initialise the total density as zeros everywhere in a 3D array of the same 
# size as the xcoord/ycoord/vcoord volumes

totalDensity <- array(data=0,dim=dim(xcoord))

# ------------------------------------------------
# Loop across moving objects

listOfObjects <- unique(dfAll$ID)
noOfObjects <- length(listOfObjects)

for (shk in 1:noOfObjects) {
  #for (shk in 1:1) {
  
  # Initialise individual object density as zeros everywhere in a 3D array of the same 
  # size as the xcoord/ycoord/vcoord volumes
  objectDensity <- array(data=0,dim=dim(xcoord))
  
  # find where this object shk is in the table
  objectIndices <- which(dfAll$ID == listOfObjects[shk])
  
  # extract only trajectories of this one object
  trajData <- dfAll[objectIndices,]
  head(trajData)
  
  # ------------------------------------------------
  # Count the number of different trajectories for this object 
  
  # listOfTrajIDs = list of unique ID values in first column of trajData 
  # listOfTrajLengths = how many points there are in each trajectory
  # trajNo = number of trajectories
  
  listOfTrajIDs <- unique(trajData$ID_traj)
  trajNo <- length(listOfTrajIDs)
  listOfTrajLengths <- vector(mode="numeric",length=trajNo)
  
  for (i in 1:trajNo) {
    # find indices where this trajectory is in dfAll table
    trajIndices <- which(trajData$ID_traj == listOfTrajIDs[i])
    
    # assign this to the vector of lengths
    listOfTrajLengths[i] <- length(trajIndices)
  }
  
  # Loop through trajData using the list of unique traj IDs
  # and at each step extract the trajectory and calculate the density.
  
  currentRow <- 1
  
  for (i in 1:trajNo) {
    # Find ID of the current trajectory
    currentTrajID <- listOfTrajIDs[i]
    currentTrajLength <- listOfTrajLengths[i]
    
    # currentLine will be a submatrix in trajData starting at currentRow with
    # currentTrajLength rows and of columns 2, 3, 4, which have coordinate
    # data
    currentLine <- trajData[currentRow:(currentRow+currentTrajLength-1),3:5]
    
    currentRow <- currentRow+currentTrajLength
    
    # Calculate density for current trajectory 
    print('*********')
    print(paste(sprintf('Object no.: %d',shk), sprintf('out of %d', noOfObjects)))
    print(paste(sprintf('Its trajectory: %d',i), sprintf('out of %d',trajNo)))
    print(paste(sprintf('Density for trajectory no.: %d',currentTrajID),sprintf('with this many points: %d', currentTrajLength) ))
    
    # Sometimes there may be only one point in a trajectory, meaning that there won't be a segment needed for
    # stacked density calculation. If that is the case, current Line will only have 1 point, so we just calculate
    # density around this one point.
    # Test data has an example of this, object 2 has trajectory 5 with only one point.
    
    if (currentTrajLength > 1) {
      
      # if the trajectory has two points or more, we calculate stacked densities around each segment
      density <- DensityAroundTrajectory(currentLine,kernelSize,xcoord,ycoord,zcoord,method,voxelSize,sigma11,sigma12) 
    
      } else {
      # if however the trajectory has only one point, then we calculate the density around this point
      density <- DensityAroundOnePoint(currentLine,kernelSize,xcoord,ycoord,zcoord,method,voxelSize,sigma11,sigma12)

          } # if (currentTrajLength > 1)
    
    # Add to object density 
    objectDensity <- objectDensity+density
    
  } # for i in trajNo
  
  # Normalise object density with no of trajectories (days)
  objectDensity  <- objectDensity/trajNo
  
  # Export density of this individual object
  objectDensity.file <- paste(c(listOfObjects[shk],"_Density.csv"),collapse="")
  WriteCSVtableFromFourRectangular3DArrays(xcoord,ycoord,zcoord,objectDensity,objectDensity.file)
  
  # Add this one individual's density to total Density
  totalDensity <- totalDensity + objectDensity
  
} # for (shk in 1:noOfObjects)

# At the end we need to normalise the total density also w/no. of objects

totalDensity  <- totalDensity/noOfObjects

# ------------------------------------------------
# Quickly test plot the results to check what happened

# Generate colour scheme
# palPurple <- brewer.pal(9,'Purples')
# colPal <- colorRampPalette(palPurple[2:9])(length(unique(totalDensity)))

# scatter3D(xcoord,ycoord,zcoord,colvar=totalDensity,pch=".",cex=2,col=colPal)

# ------------------------------------------------
# Export density into csv file that can be imported into Voxler.

print('Exporting total density into a csv file.')
WriteCSVtableFromFourRectangular3DArrays(xcoord,ycoord,zcoord,totalDensity,output.file)

# ------------------------------------------------


