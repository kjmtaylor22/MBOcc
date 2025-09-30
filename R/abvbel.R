#' @title Recode community table to above (1) or below (0) species mean
#'
#' @description These functions facilitate data carpentry for community table and metadata.
#' @param x To be used with `apply`. Single element of community table
#' @param stat The statistic to be calculated.
#' @keywords above-or-below-mean
#' @export
#' @examples
#' comm <- apply(comm, 2, u.d)



u.d <- function(x, stat=c("mean","logmean","median")){ ## function to convert table to 1/0 for above/below colMeans, respectively
  if (stat=="mean"){
    y <- mean(x) #exclude NAs from mean calculation
    test <- expression(x1<=y)
  }
  if (stat=="logmean"){
    y <- mean(x, na.rm=T) #exclude NAs from mean calculation
    x[is.na(x)] <- 0 #return NAs to 0s
    test <- expression(x1<y)
  }
  if (stat=="median"){
    y <- median(x) #exclude NAs from mean calculation
    test <- expression(x1<=y)
  }
  z <- sapply(x, FUN=function(x1){if(eval(test)|is.nan(y)){return(0)}else{return(1)}})
  return(z)
}
