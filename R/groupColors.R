#' @title Group color palette by taxonomy
#' @description Generate a rainbow color palette for taxa, grouped by genus or lower.
#' @author Kara J.M. Taylor (`k.taylor2@ufl.edu`)
#' @param taxa Character string of taxon names as given in MBOcc.obj
#' @export



group.colors <- function(taxa, show=T){
  library(grDevices)
  ramp <- c("magenta", "red", "orange", "yellow", "green",
            "lightseagreen", "mediumblue", "darkorchid4")
  a <- strsplit(taxa, split = " ")
  b <- unlist(lapply(a, FUN=function(x){return(x[[1]])}))
  c <- summary(as.factor(b))
  clr <- colorRampPalette(ramp)(length(c))
  out <- c()
  for (i in 1:length(c)){
    d <- colorRampPalette(c(clr[i],"white"))(c[i]+2)
    out <- c(out, d[1:c[i]])
  }
  names(out) <- sort(taxa)
  if (show==T){scales::show_col(out)}
  return(out)
}
