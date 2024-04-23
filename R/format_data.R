#' @title Format multi-species data for multi-site occupancy
#' @description
#' A short description...
#' @author Jose Miguel Ponciano Castellanos (`josemi@ufl.edu`), Kara J.M. Taylor (`k.taylor2@ufl.edu`)
#' @param
#' @param
#' @param
#' @export

format <- function(comm, meta, tax, zeroes){

  library(dplyr)

  if ("datatables"%in%class(tax)){
    taxonomy <- function(taxrow){
      suppressWarnings(i <- min(which(is.na(taxrow)))-1)
      if (is.infinite(i)){
        i <- paste(taxrow[length(taxrow)-1], taxrow[length(taxrow)])
        return(i)
      } else {
        return(taxrow[i])
      }
    }
    tax <- tax$x$data
    tax <- cbind(tag=tax[,1], tax[,-1], taxonomy=apply(tax, 1, taxonomy)) %>%
      mutate(otu.name=paste0(taxonomy, " (", tag, ")"))
  }

  if (nrow(comm)!=nrow(meta)|ncol(comm)!=nrow(tax)){
    stop("Check dataframe dimensions for congruence.")
  }

  if (!any(row.names(comm)%in%row.names(meta))){
    stop("Check row names of community and metadata matrices for identity.")
  }

  if (!any(colnames(comm)%in%row.names(tax))){
    stop("Check community column names and taxonomy row names identity.")
  }

  # change community table to presence-absence
  comm[comm>0] <- 1

  # make sure community matrix order matches taxonomy and metadata
  comm <- comm[row.names(meta),]
  comm <- comm[,row.names(tax)]





}
