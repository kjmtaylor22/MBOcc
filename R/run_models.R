#' @title Apply occupancy model to formatted MB list
#' @description
#' A short description...
#' @author Jose Miguel Ponciano Castellanos (josemi@ufl.edu), Kara J.M. Taylor (`k.taylor2@ufl.edu`)
#' @param MBOcc.obj Formatted object output by `format` function.
#' @param formulae Vector of formulae to run for each psi.
#' Should be of length equal to `assign.psi`.
#' @param assign.psi List of numeric strings indicating how many psi to estimate.
#' List should be of length equal to `formulae`.
#' @export

MBOcc <- function(MBOcc.obj, formulae, assign.psi){

  MBruntime(n.sp=length(MBOcc.obj),
            n.form=length(formulae),
            n.psi=length(assign.psi))

  library(dplyr)

  get.groups <- function(x){
    try <- grep("\\b\\d{4}\\b", names(x), invert = T, value=T)
    return(try)
  }

  run.models <- function(x){
    groups <- get.groups(x)
    out <- vector("list", length(assign.psi))
    for (i in 1:length(assign.psi)){
      for (j in 1:length(formulae)){
        try <- formulae[[j]]
        for (f in 2:length(assign.psi[[i]])){try <- append(try, formulae[[j]])}
        if (try[[1]]=="~1"){
          out[[i]][[j]] <- const.model(formula=try, data=x, assign.psi=assign.psi[[i]])
        } else {
          for (g in groups){
            if (is.numeric(x[,g])){
              if (any(grepl("*", try[[1]][[2]], fixed=T))){next}
              out[[i]][[j]] <- const.model(formula=try, data=x, assign.psi=assign.psi[[i]])
            } else {
              if(any(grepl(g, try[[1]][[2]]))){
                check <- length(unique(x[,g]))
                if (check>=2){
                  out[[i]][[j]] <- const.model(formula=try, data=x, assign.psi=assign.psi[[i]])
                }
              }
            }
          }
        }
      }
    }
    return(out)
  }

  min.BIC <- function(x, second.filter=F){
    if (any(unlist(lapply(x, is.null)))){
      x <- x[-which(unlist(lapply(x, is.null)))]}
    y <- as.numeric(unlist(x)[grep("BIC", names(unlist(x)))])
    miny <- min(y)
    deltaBIC <- sort(y-miny)[2]
    z <- x[[which(y==miny)]]
    if (second.filter==T){
      deltaBIC <- min(z$DeltaNextBestModel, deltaBIC)}
    if (length(y)>1){z$DeltaNextBestModel <- deltaBIC}
    return(z)
  }

  tabulate.models <- function(x){
    call <- paste(x$call, collapse=", ")
    psi <- paste(x$assigned.psi, collapse=".")
    bic <- x$BIC
    deltaNext <- x$DeltaNextBestModel
    out <- data.frame(formulae=call, BIC=bic, psi=psi, deltaNext=deltaNext)
    return(out)
  }



  ## run models and format output
  model.list2 <- lapply(MBOcc.obj, run.models)

  models.out <- list()
  for(i in names(MBOcc.obj)){
    models.out[[i]] <- lapply(model.list2[[i]], min.BIC)
  }
  models.out2 <- lapply(models.out, min.BIC, second.filter=T)

  models.table <- data.frame()
  tmp <- lapply(models.out2, FUN=tabulate.models)

  m.tmp <- data.frame(stringsAsFactors = F)
  for (j in 1:length(tmp)){
    m.tmp <- rbind(m.tmp, tmp[[j]])
  }

  models.table <- data.frame(taxon=names(models.out2), m.tmp)

  outbag <- list(estimates=models.out2, best=models.table)
  class(outbag) <- "MBOcc.obj.est"
  return(outbag)
}
