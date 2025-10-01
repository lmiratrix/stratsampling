
# This allows passing either a vector of values or a character name
# for a column for a more flexible command
get_col_or_vector = function(dat, input) {
  name <- deparse(substitute(input))
  if ( is.null(dat) ) {
    stopifnot( is.vector(input) )
    return( input )
  } else if (length(input) == 1) {
    if ( is.character(input) ) {
      if ( !(input %in% colnames(dat)) ) {
        stop(paste("Column", input, "not found in data frame"))
      }
    }
    return(dat[[input]])
  } else if (length(input) == nrow(dat)) {
    return(input)
  } else {
    stop(paste("Input vector ", name, " must be either a column name or a vector of length ", nrow(dat)) )
  }
}
