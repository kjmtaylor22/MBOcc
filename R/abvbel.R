#' @title Recode community table to above (1) or below (0) species mean
#'
#' @description These functions facilitate data carpentry for community table and metadata.
#' @param x To be used with `apply`. Single element of community table
#' @keywords above-or-below-mean
#' @export
#' @examples
#' comm <- apply(comm, 2, u.d)



u.d <- function(x){ ## function to convert table to 1/0 for above/below colMeans, respectively
  y <- mean(x)
  z <- sapply(x, FUN=function(x1){if(x1>y){return(1)}else{return(0)}})
  return(z)
}
