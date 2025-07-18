# a function to convert bits to integers/numeric
#' *NOTE: this function uses R-type (i.e., reversed) bits!*
bit_to_int <- function(bit) {
  if(class(bit) == 'raw') {
    bit <- as.numeric(bit)
  } else if(class(bit) == 'character') {
    bit <- as.numeric(strsplit(bit, '*')[[1]])
    } else if(! class(bit) %in% c('numeric', 'integer')) {
    stop('`bit` must be of class "raw", "character", numeric" or "integer".\n')
  }
  
  if(length(bit) != 16) {
    stop('This function can only deal with 16-bit signed shorts.')
  }
  
  as.numeric(t(2^(0:14)) %*% bit[1:15] - 2^15 * bit[16])
}

if(FALSE) { # some tests
  bit_to_int(intToBits(-10))
  bit_to_int(intToBits(0))
  bit_to_int(intToBits(10))
}