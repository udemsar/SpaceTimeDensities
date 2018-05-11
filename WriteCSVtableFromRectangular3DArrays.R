# ------------------------------------------------
# WriteCSVtableFromFourRectangular3DArrays.R
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

# ------------------------------------------------
# Supporting file to StackedSpaceTimeDensities.R
# ------------------------------------------------
# Writes data from the four 3D arrays representing xcoordinates, ycoordinates,
# zcoordinates and vcoordinates into a CSV file which can
# be imported to Voxler. CSV has the following four columns: xcoord,
# ycoord, zcoord (time in seconds) and vcoord (density value).
# ------------------------------------------------

WriteCSVtableFromFourRectangular3DArrays <- function(xc,yc,zc,vc,FileName) {
  
  dimensions <- dim(xc)
  rowsX <- dimensions[1]
  colsX <- dimensions[2]
  sheetsX <- dimensions[3]
  
  # Build an empty CSV table of correct dimensions
  CSVtable <- array(data=0,dim=c(rowsX*colsX*sheetsX,4))
  colnames(CSVtable) <- c('X','Y','T','density')
  
  # Convert the volumetric representation into x,y,t,density table
  for (k in 1:sheetsX) {
    #print('.');
    for (i in 1:rowsX) { 
      for (j in 1:colsX) {
      # print(c(i,j,k))
      CSVtable[(k-1)*rowsX*colsX+(i-1)*colsX+j,] <- c(xc[i,j,k],yc[i,j,k],zc[i,j,k],vc[i,j,k])
      } # j
    } # i
  } # k 
  
  # export into the given output file
  write.csv(CSVtable,FileName,row.names=FALSE)
  
  #return(f)  
}