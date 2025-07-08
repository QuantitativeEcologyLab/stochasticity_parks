# constant to move values away from boundaries
ndvi_to_01 <- function(ndvi, avoid_boundaries = FALSE) {
  ndvi_scaled <- (ndvi + 1) / 2
  
  if(avoid_boundaries) {
    EPS <- .Machine$double.eps * 100
    
    ndvi_scaled <- dplyr::if_else(condition = ndvi_scaled - EPS <= 0,
                                  true = EPS,
                                  false = ndvi_scaled)
    ndvi_scaled <- dplyr::if_else(condition = ndvi_scaled + EPS >= 1,
                                  true = ndvi_scaled - EPS,
                                  false = ndvi_scaled)
  }
  
  return(ndvi_scaled)
}

ndvi_to_11 <- function(ndvi_scaled, avoided_boundaries = FALSE) {
  if(avoided_boundaries) {
    EPS <- .Machine$double.eps * 100
    
    ndvi_scaled <- dplyr::if_else(condition = ndvi_scaled == EPS,
                                  true = 0, # = ndvi_scaled - EPS,
                                  false = ndvi_scaled)
    ndvi_scaled <- dplyr::if_else(condition = ndvi_scaled == 1 - EPS,
                                  true = 1, # ndvi_scaled + EPS,
                                  false = ndvi_scaled)
  }
  
  ndvi <- ndvi_scaled * 2 - 1
  
  return(ndvi)
}
