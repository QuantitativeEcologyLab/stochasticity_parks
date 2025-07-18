#' to drop values flagged with quality flag assurance `.qa`. Note that the
#' flag positions are counted starting from 0 for consistency. For example,
#' if `flag = 2` and we want to exclude values (10, 11) and keep (00, 11),
#' we can use `is_flagged() == (QA %% (2 * 2) >= 2)`:
#' 00 keep; return TRUE  to drop
#' 01 keep; return TRUE  to drop
#' 10 drop; return FALSE to keep
#' 10 drop; return FALSE to keep
#' if you want to do the opposite, just negate the output with `!`
#' `QA` can be a single value or a SpatRaster
source('functions/bit_to_int.R')

is_flagged <- function(QA, flag_position) {
  if(length(QA) > 1) {
    stop('`QA` must be of length 1.\n')
  }
  if(length(flag_position) > 1) {
    stop('`flag_position` must be of length 1.\n')
  }
  
  # find bit value using powers (incices start from 0)
  return(QA %% ((2^flag_position) * 2) >= (2^flag_position))
}
