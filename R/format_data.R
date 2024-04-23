#' @title Format multi-species data for multi-site occupancy
#' @description
#' A short description...
#' @author Jose Miguel Ponciano Castellanos (`josemi@ufl.edu`), Kara J.M. Taylor (`k.taylor2@ufl.edu`)
#' @param
#' @param
#' @param
#' @param zeroes Minimum number of samples from which a species can be absent.
#' @export

format <- function(comm, meta, tax, zeroes, states){

  library(dplyr)

  reframe <- function(x, states){
    y <- dplyr::mutate(x, state=paste(contains(unique(states)), sep=""))
    z <- reshape2::dcast(y, . ~ state)
    return(z)
  }

  add.states <- function(x){
    states <- names(trans.list[[1]])
    missing.states <- which(!states%in%names(x))
    add.states <- matrix(0,nrow(x),length(missing.states),dimnames=list(NULL,states[missing.states]))
    y <- cbind(x, add.states)
    return(y)
  }



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
      mutate(tag.name=paste0(taxonomy, " (", tag, ")"))
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

  # remove species for which they are absent more than the threshold in param `zeroes`
  pull <- apply(comm, 2, function(x){length(x[which(x==0)])})
  comm <- comm[,-which(pull>zeroes)]

  # reduce taxonomy to match reduced community table
  tax <- tax[colnames(comm),]

  # create a new data frame with the reduced information
  taxa <- t(comm) %>% as.data.frame() %>%
    tibble::rownames_to_column(var="tag") %>%
    left_join(tax[,c("tag","taxonomy","tag.name")]) %>%
    as.data.frame() %>% `row.names<-`(.$tag.name) %>%
    .[,-c(1,(ncol(.)-2):ncol(.))] %>% t()

  # make sure SampleID is specified in the metadata
  meta <- data.frame(SampleID=row.names(meta), meta)

  comm.list <- list()
  for (i in colnames(taxa)){

    tax.comm <- data.frame(SampleID=row.names(taxa), Taxon=taxa[,i]) %>%
      left_join(meta) %>% droplevels(.)

    if (!is.null(states)){
      tax.comm <- reshape2::dcast(tax.comm, . ~ vars(states), value.var = "Taxon")
    }

    if (any(is.na(tax.comm))){tax.comm[is.na(tax.comm)] <- 0}

    comm.list[[i]] <- tax.comm
  }

  if (!is.null(states)){
    trans.list <- lapply(comm.list, reframe, states=states)
    stand.list <- lapply(trans.list, add.states)
  }
  return(stand.list)
}
