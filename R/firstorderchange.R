#' @title Recode timeseries community table to first order sequential change
#'
#' @description These functions facilitate data carpentry for community table and metadata.
#' Based on the abundance at time `t`, recode the abundance at time `t+1`:
#' `N_{t+1} > N_{t} = 1` or `N_{t+1} < N_{t} = 0`
#' @param comm Community matrix
#' @param id ID variable indicating unit of repeated measurements
#' @param sitevar vector indicating which samples represent time-repeated measures
#' @param timevar vector indicating the time step at which the samples were taken
#' @param sampleID vector indicating the sample IDs
#' @keywords first-order sequence
#' @export
#' @examples
#' comm <- apply(comm, 2, f.o, sitevar=meta$site, timevar=meta$time)



f.o <- function(comm, id, sitevar, timevar, sampleID){
  library(tibble)
  library(dplyr)
  library(reshape2)

  suppressMessages(

  comm <- lapply(comm, function(z){

    df <- data.frame(id=id, sitevar=sitevar, timevar=timevar, z=z, row.names=sampleID)
    df$id <- as.character(df$id)
    if (!is.numeric(df$timevar)){
      if (class(df$timevar)=="factor"){df$timevar <- as.character(df$timevar)}
      df$timevar <- as.numeric(df$timevar)
    }

    df <- df[order(df$id, df$sitevar, df$timevar),]

    cast <- acast(df, timevar ~ sitevar ~ id, value.var="z")
    cast <- apply(cast, 3, identity, simplify = FALSE)
    cast <- lapply(cast, na.omit)

    to.fo <- function(x){
      binary <- matrix(0,nrow=nrow(x)-1, ncol=ncol(x))
      for (t in 2:nrow(x)){
        binary[t-1,] <- x[t,]-x[t-1,]
      }
      binary[binary<0] <- 0
      binary[binary>0] <- 1

      remelt <- as.data.frame(binary) %>%
        `colnames<-`(colnames(x)) %>%
        data.frame(timevar=row.names(x)[-1]) %>%
        melt(variable.name="sitevar", id.vars="timevar")

      return(remelt)
    }

    out <- lapply(cast, to.fo) %>% melt()
    out$L1 <- as.character(out$L1)

    df$timevar <- as.character(df$timevar)
    df$SampleID <- row.names(df)

    df <- left_join(out, df, by=c("timevar"="timevar","sitevar"="sitevar","L1"="id"))
    out <- df[,c(4, ncol(df))]
    return(out)

  })

  )

  suppressMessages(
  out <- reshape2::melt(comm) %>%
    reshape2::dcast(SampleID ~ L1) %>%
    column_to_rownames("SampleID")
  )

  return(out)
}
