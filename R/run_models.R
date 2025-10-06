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

MrSpOcc <- function(MBOcc.obj, formulae, assign.psi){

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
    out <- vector("list", length(formulae))
    for (i in 1:length(assign.psi)){
      for (j in 1:length(formulae)){
        try <- formulae[[j]]
        for (f in 2:length(assign.psi[[i]])){try <- append(try, formulae[[j]])}
        if (try[[1]]=="~1"){
          out[[j]][[i]] <- const.model(formula=try, data=x, assign.psi=assign.psi[[i]])
        } else {
          for (g in groups){
            if (is.numeric(x[,g])){
              if (any(grepl("*", try[[1]][[2]], fixed=T))){next}
              out[[j]][[i]] <- const.model(formula=try, data=x, assign.psi=assign.psi[[i]])
            } else {
              if(any(grepl(g, try[[1]][[2]]))){
                check <- length(unique(x[,g]))
                if (check>=2){
                  out[[j]][[i]] <- const.model(formula=try, data=x, assign.psi=assign.psi[[i]])
                }
              }
            }
          }
        }
      }
    }
    return(out)
  }

  min.BIC <- function(x){
    if (any(unlist(lapply(x, is.null)))){
      x <- x[-which(unlist(lapply(x, is.null)))]}

    y <- as.numeric(unlist(x)[grep("BIC", names(unlist(x)))])
    miny <- min(y)
    diffs <- y-miny

    deltaBIC <- sort(diffs[-which(diffs==0)])[1]
    deltaNext <- paste(paste(x[[which(diffs==deltaBIC)[1]]]$assigned.psi, collapse="."),
                       x[[which(diffs==deltaBIC)[1]]]$call[1], sep=": ")

    z <- x[[which(y==miny)]]
    if (length(y)>1){
      z$DeltaNextBIC <- deltaBIC
      z$DeltaNextBestModel <- deltaNext
    }

    return(z)
  }

  tabulate.models <- function(x){
    call <- paste(x$call, collapse=", ")
    psi <- paste(x$assigned.psi, collapse=".")
    bic <- x$BIC
    deltaBIC <- x$DeltaNextBIC
    deltaNext <- x$DeltaNextBestModel
    out <- data.frame(formulae=call, BIC=bic, psi=psi,
                      deltaBIC=deltaBIC, deltaNext=deltaNext)
    return(out)
  }

  ## run models and format output
  model.list2 <- lapply(MBOcc.obj, run.models)

  models.out <- lapply(model.list2, unlist, recursive=F)

  for (i in 1:length(models.out)){try <- min.BIC(models.out[[i]])}

  models.out2 <- lapply(models.out, min.BIC)

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
