# -0.1 is the default lower boundary bc it's the lowest value in the data
ndvi_to_01 <- function(ndvi, lower = -0.1, upper = 1) {
  ndvi_scaled <- (ndvi - lower) / (upper - lower)
  
  return(ndvi_scaled)
}

ndvi_to_11 <- function(ndvi_scaled, lower = -0.1, upper = 1) {
  ndvi <- ndvi_scaled * (upper - lower) + lower
  
  return(ndvi)
}

if(FALSE) {
  # test the functions
  ndvi_to_01(c(-0.1, 0, 1))
  ndvi_to_11(c(0, 1))
  ndvi_to_11(ndvi_to_01(c(-0.1, 0, 1)))
}
