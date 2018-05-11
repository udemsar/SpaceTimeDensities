# ------------------------------------------------
# Decay functions for the calculations of space-time density
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
# This file includes three decay functions:
# - linear
# - bisquare
# - Gaussian
# ------------------------------------------------


# ------------------------------------------------
# 1. Linear decay function
# ------------------------------------------------
# Calculates the linear decay weighting function:
# f(x): [0,h] -> [0,1], so that 0->1 and h->0
# f(x)=(h-x)/h
# h must be positive (kernel distance), x must be in [0,h]
#
# Function is scaled as 2D probability, so that the integral
# under decay function over each voxel layer is equal to 1.
# f=(3/(pi*h^2))*(h-x)/h;
#
# To get the density value, this is then multiplied with voxel area.
# ------------------------------------------------

LinearDecay <- function(x,h,vxl) {
  
# Calculate linear decay in 2D
f <- (3/(pi*h^2))*(h-x)/h

# Scale to voxel area:
# multiply with the 2D surface of the voxel
f <- f*vxl^2

return(f)
}

# ------------------------------------------------
# Bisquare decay function
# ------------------------------------------------
# ------------------------------------------------
# Calculates the appropriate bisquare weighting function:
# f(x): [-h,h] -> [0,1], where -h and h go to 0 and 0 goes to 1.
# f(x)=1/h^4 * (x^2-h^2)^2
# h must be positive (kernel distance), x must be in [-h,h]
#
# Function is scaled as 2D probability, so that so that the integral
# under decay function over each voxel layer is equal to 1.
# f=(3/(pi*h^2))*(h-x)/h;
#
# To get proper density, this is then multiplied with voxel area.
# ------------------------------------------------

BisquareDecay <- function(x,h,vxl) {
  
  # Bisquare decay in 2D
  f <- (3/(pi*h^2))*(1/h^4) * (x^2-h^2)^2
  
  # Scale to voxel area:
  # multiply with the 2D surface of the voxel
  f <- f*vxl^2
  
  return(f)
}


# ------------------------------------------------
# Gaussian decay function
# ------------------------------------------------
# Calculates the appropriate Gaussian decay weighting function (bell):
# f(x)=1/sqrt(2*pi*sigma2) * exp(-(x-mu)^2/(2*sigma2))
# mu = mean, location of the peak of the function
# sigma2 = variance, width of the bell
# 
# Function is scaled as 2D probability, so that so that the integral
# under decay function over each voxel layer is equal to 1.
# f=1/(2*pi*sigma2) * exp(-(x-mu)^2/(2*sigma2));
#
# To get proper density, this is then multiplied with voxel area.
# ------------------------------------------------

GaussianDecay <- function(x,mu,sigma2,vxl) {

# Gaussian decay in 2D
f=1/(2*pi*sigma2) * exp(-(x-mu)^2/(2*sigma2));

# Scale to voxel area:
# multiply with the 2D surface of the voxel
f <- f*vxl^2

return(f)
}