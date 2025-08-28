#' @title Combine lists with the same names to append into new sub-lists
#' @description Match lists by name and append
#' @author Kara J.M. Taylor (`k.taylor2@ufl.edu`)
#' @param list1 List to match to
#' @param list2 List to match
#' @param FUN Appending function, e.g. `append`, `c`, `rbind`, `cbind`
#' @param recursive Logical; Should list elements also be checked for concatenation?
#' @export


recombine <- function(list1, list2, FUN, recursive=F){

  FUN <- match.fun(FUN)

  #copy list1 to output so it becomes the base
  output <- list1
  #identify names present in both lists
  common_names <- intersect(names(list1), names(list2))
  #rbind elements with common names
  output[common_names] <- Map(FUN, list1[common_names], list2[common_names])
  #get unique names from list2
  unique_names_list2 <- setdiff(names(list2), names(list1))
  #Add unique elements from list2
  output[unique_names_list2] <- list2[unique_names_list2]

  if (recursive){
    output <- lapply(output, function(x){
      detect <- any(duplicated(names(x)))
      if (detect){
        pull <- which(duplicated(names(x)))
        y <- x[pull]; x <- x[-pull]
        match <- match(names(y), names(x))
        x[match] <- Map(FUN, x[match], list(y))
      }
      return(x)
    })
  }
  return(output)
}
