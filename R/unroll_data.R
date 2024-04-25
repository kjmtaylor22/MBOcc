#' @title Unroll MBOcc estimates from model output
#' @description Format MBOcc output into table format for plotting
#' @author Kara J.M. Taylor (`k.taylor2@ufl.edu`)
#' @param
#' @export

unroll <- function(MBOcc.obj, MBOcc.est, group.vars, site.list){
  states <- names(MBOcc.obj[[1]])[!names(MBOcc.obj[[1]])%in%group.vars]

  outlist <- lapply(MBOcc.est$estimates, function(x){
    x <- x[c(7,9)];
    colnames(x[[1]]) <- states;
    colnames(x[[2]]) <- site.list;
    y <- lapply(x,function(z){
      return(cbind(MBOcc.obj[[1]][,group.vars], z))})})

  outtable <- reshape2::melt(outlist, id.vars=group.vars)
  outtable <- dplyr::left_join(outtable, MBOcc.est$best, by=c("L1"="taxon"))
  outtable$L1 <- as.factor(outtable$L1)

  outbag <- list(P.site=subset(outtable, subset=L2=="regs.mat"),
                 P.state=subset(outtable, subset=L2=="preds.mat"))

  class(outbag) <- "MBOcc"
  return(outbag)
}
