#' @title Combine lists with the same names to append into new sub-lists
#' @description Match lists by name and append
#' @author Kara J.M. Taylor (`k.taylor2@ufl.edu`)
#' @param list1 List to match to
#' @param list2 List to match
#' @param FUN Appending function, e.g. `append`, `c`, `rbind`, `cbind`
#' @param recursive Logical; combine lower list elements of single list
#' @export


recombine <- function(list1, list2=NULL, FUN, recursive=F){

  FUN <- match.fun(FUN)

  doMap <- function(X,Y){
    #copy list1 to output so it becomes the base
    output <- X
    #identify names present in both lists
    common_names <- intersect(names(X), names(Y))
    #rbind elements with common names
    output[common_names] <- Map(FUN, X[common_names], Y[common_names])
    #get unique names from list2
    unique_names_list2 <- setdiff(names(Y), names(X))
    #Add unique elements from list2
    output[unique_names_list2] <- Y[unique_names_list2]

    return(output)
  }

  if (!recursive){

    stopifnot(!is.null(list2))
    out <- doMap(list1, list2)

  } else {

    check <- unlist(lapply(list1, is.list))
    stopifnot(length(check)>1&all(check))

    n <- length(list1)
    out <- list1[[1]]
    for (i in 2:n){
      out <- doMap(out, list1[[i]])

    }
  }
  return(out)
}
