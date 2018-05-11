# ------------------------------------------------
# SquaredDistance2points.m
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
# Calculates squared distance between 2 points in 3D, 
# given as (x,y,z) and (x1,y1,z1).
# ------------------------------------------------

SquaredDistance2points <- function(x,y,z,x1,y1,z1) {
  
  f <- (x-x1)^2+(y-y1)^2+(z-z1)^2

  return(f)  
  }