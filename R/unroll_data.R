#' @title Unroll MBOcc estimates from model output
#' @description Format MBOcc output into table format for plotting
#' @author Kara J.M. Taylor (`k.taylor2@ufl.edu`)
#' @param MBOcc.obj Formatted MBOcc input data list
#' @param MBOcc.est MBOcc output model estimates and best fit data
#' @param group.vars Character string of names of all covariates used in models
#' @param site.list Character string of site names, in the same order as given in `assigned.psi`
#' @export

unroll <- function(MBOcc.obj, MBOcc.est, group.vars, site.list){

  stopifnot(class(MBOcc.obj)=="MBOcc.obj", class(MBOcc.est)=="MBOcc.obj.est")

  library(dplyr)

  findState <- function(x){
    x <- names(x)
    sp <- strsplit(x, "")
    suppressWarnings(nu <- lapply(sp, as.numeric))
    nu <- unlist(lapply(nu, function(z){any(is.na(z))}))
    x[which(nu==F)]
  }

  states <- findState(MBOcc.obj[[1]])

  outlist <- lapply(MBOcc.est$estimates, function(x){
    x <- x[c(7,9)]
    colnames(x[[1]]) <- states
    colnames(x[[2]]) <- site.list
    y <- lapply(x,function(z){
      return(cbind(MBOcc.obj[[1]][,group.vars], z))})})

  outlist2 <- lapply(MBOcc.est$estimates, function(x){
    m <- x[c(1,2,4:6,10)]
    h <- m[[5]]
    if (nrow(h)>1){
      library(MASS)
      ih <- ginv(h)
      dih <- diag(ih)
      sdih <- sqrt(dih)
    } else {sdih <- sqrt(1/h)}
      out <- data.frame(formulae=paste(m[[1]],collapse=" "),
                        psi=paste(m[[2]], collapse=" "),
                        BIC=m[[4]], dBICnext=m[[6]],
                        mle=m[[3]], se=sdih)
      return(out)})

  outtable <- reshape2::melt(outlist, id.vars=group.vars)
  outtable <- dplyr::left_join(outtable, MBOcc.est$best, by=c("L1"="taxon"))
  outtable$L1 <- as.factor(outtable$L1)

  outlist2 <- lapply(outlist2, function(x){
    f <- strsplit(x$formulae, "~")[[1]][-1]
    f <- paste("~",f)
    f <- gsub(" ", "", f)
    p <- strsplit(x$psi, " ")[[1]]
    p1 <- unique(cbind(f, p))
    betas <- data.frame()
    for (i in 1:nrow(p1)){
      mm <- colnames(model.matrix(as.formula(p1[i,1]), data=MBOcc.obj[[1]]))
      suppressWarnings(betas <- rbind(betas, data.frame(site.psi=p1[i,2], site.call=p1[i,1], beta=mm)))
    }
    out <- cbind(betas, model.psi=x$psi, model.call=x$formulae, mle=x$mle, se=x$se)
    return(out)
  })

  sigstar <- function(x){
    out <- ""
    if (x < 0.1){out <- "."}
    if (x < 0.05){out <- "*"}
    if (x < 0.01){out <- "**"}
    if (x < 0.001){out <- "***"}
    return(out)
  }

  outtable2 <- suppressMessages(reshape2::melt(outlist2)) %>%
    unique() %>%
    reshape2::dcast(model.psi+model.call+site.psi+site.call+L1+beta ~ variable) %>%
    mutate(ucl=mle+se*1.96, lcl=mle-se*1.96, statT=(mle-0)/se) %>%
    mutate(p=signif(2*pnorm(abs(statT), lower.tail = FALSE), digits = 3)) %>%
    cbind(sig=vapply(.$p, sigstar, "vector"))


  outbag <- list(P.site=droplevels(subset(outtable, subset=L2=="regs.mat")),
                 P.state=droplevels(subset(outtable, subset=L2=="preds.mat")),
                 MLE.SE=outtable2,
                 groups=group.vars,
                 sites=site.list)

  class(outbag) <- "MBOcc"
  return(outbag)
}
