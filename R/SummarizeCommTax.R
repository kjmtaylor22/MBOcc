#' @title Summarize count data by higher taxon
#'
#' @description Adds up counts of taxon columns by some higher taxon identification. The second data frame should provide a key for mapping the column IDs to higher taxon names.
#' @param comm community matrix
#' @param tax full matrix containing 7 taxonomic levels (out from SILVA reference)
#' @param level column of traits matrix containing taxonomic name to summarize to
#' @export
#'
#'

SummarizeCommTax <- function(comm, tax, level=NULL){

  library(dplyr)

  flip <- t(comm) %>% as.data.frame()
  flip <- data.frame(tag=row.names(flip), flip)

  if (!is.null(level)){

    collapse1 <- left_join(flip, tax[,c("tag", "taxonomy", level)]) %>% .[,-1] %>%
      group_by_("taxonomy", level) %>% summarize_all(funs(sum)) %>% as.data.frame()

    collapse2 <- collapse1[collapse1[,2]!="",] %>% .[,-1] %>%
      group_by_(level) %>% summarize_all(funs(sum)) %>% as.data.frame()

    collapse1 <- collapse1[collapse1[,2]=="",] %>% .[,-2] %>%
      mutate(taxonomy=paste("unclassified", taxonomy))

    colnames(collapse2)[1] <- "taxonomy"

    collapse <- rbind(collapse2, collapse1)

  } else {

    collapse <- left_join(flip, tax[,c("tag", "taxonomy")]) %>% .[,-1] %>%
      group_by(taxonomy) %>% summarize_all(funs(sum)) %>% as.data.frame()

  }

  row.names(collapse) <- collapse$taxonomy

  collapse <- collapse[,-1] %>% t() %>% as.data.frame()

  return(collapse)
}

TaxSummarize <- function(comm, traits){

  library(dplyr)

  flip <- t(comm) %>% as.data.frame()
  flip <- data.frame(tag=row.names(flip), flip)

  if (any(is.na(traits))){
    n <- which(is.na(traits))
    traits[n] <- "unidentified"
  }

  collapse <- data.frame(tax=traits, flip) %>%
    .[,-2] %>% group_by(tax) %>% summarize_all(funs(sum)) %>%
    as.data.frame() %>% `row.names<-`(.$tax) %>% .[,-1] %>%
    t() %>% as.data.frame()

  return(collapse)
}


