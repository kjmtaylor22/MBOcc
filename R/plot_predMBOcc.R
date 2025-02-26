#' @title Plot predictions from predict.MBOcc.obj.est
#' @description Visualize predictions
#' @author Kara J.M. Taylor (`k.taylor2@ufl.edu`)
#' @param predMBOcc Output list from `pred.MBOcc` containing newdata and predicted values
#' @param covs List of covariates in the output models
#' @param facets List of faceting variables in the output models
#' @param return Boolean; return all `ggplot` as list
#' @export

plot.predMBOcc <- function(predMBOcc, covs, facets, return=T){
  library(ggplot2)
  library(egg)
  library(ggpubr)

  stopifnot(class(predMBOcc)=="predMBOcc")

  colsSP <- group.colors(names(predMBOcc), show=F)

  out <- list()

  if (return==F){
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  #try <- lapply(predMBOcc, function(x){if (!"Sites"%in%names(x)){cbind(x,Sites=rep(1,nrow(x)))}else{return(x)}})
  try <- predMBOcc
  class(try) <- "list"

  out <- list()

  if (!is.null(facets)){
    columns <- paste(paste(c("formulae","psi","L1","Sites",covs,facets), collapse="+"), "~ variable")
    preds <- reshape2::melt(try, id.vars=c("formulae","psi","Sites",covs,facets)) %>%
      unique() %>%
      reshape2::dcast(as.formula(columns))
  } else {
    columns <- paste(paste(c("formulae","psi","L1","Sites",covs), collapse="+"), "~ variable")
    preds <- reshape2::melt(try, id.vars=c("formulae","psi","Sites",covs)) %>%
      unique() %>%
      reshape2::dcast(as.formula(columns))
  }

  check <- function(x, against){
    if (length(grep(x,against))!=0){return(T)} else {return(F)}
  }


  for (j in unique(preds$formulae)){
    if (j == "~1 ~1 ~1 ~1"){next}
    out[[j]] <- list()
    sub <- subset(preds, subset=formulae==j)

    pull <- unlist(lapply(covs,j,FUN=check))
    if (length(which(pull))!=0){
      i <- covs[which(pull)]
    } else {next}

    #for (i in covs){
    #  pull <- grep(i, j)
    #  if (length(pull)!=0){break}
    #}

    if (!is.null(facets)){
      for (f in facets){
        pull <- grep(f, j)
        if (length(pull)!=0){break}
      }
    }

    for (k in unique(sub$psi)){

      if (!is.null(facets)&length(pull)!=0){
        fm <- as.formula(paste("Sites ~", f))
        ffs <- expression(facet_grid(fm))
      } else {ffs <- expression(facet_grid(Sites~.))}

      sub2 <- subset(sub, subset=psi==k)

      for (l in i){
        g <- ggplot(sub2) +
          labs(subtitle=k, x=l, color=NULL, fill=NULL,
               y=expression("Predicted change in " ~psi)) +
          eval(ffs) +
          geom_ribbon(aes(x=eval(parse(text=l)), ymin=lcl,
                          ymax=ucl, fill=L1), alpha=0.3) +
          geom_line(aes(x=eval(parse(text=l)), y=fit, color=L1), linewidth=1) +
          scale_color_manual(values=colsSP) +
          scale_fill_manual(values=colsSP) +
          coord_cartesian(ylim=c(0,1)) +
          ggthemes::theme_few() +
          theme(legend.position = "bottom",
                legend.direction = "vertical",
                strip.text = element_text(size=7),
                axis.text = element_text(size=7))

        out[[j]][[paste(k,l,sep=":")]] <- g
      }
    }
    gg <- egg::ggarrange(plots=out[[j]], nrow=1, draw=F)
    ggg <- annotate_figure(gg, top=text_grob(j, face = "bold", size = 10))
    if (return==T){out[[j]] <- ggg} else {print(ggg)}
  }

  if (return==T){return(out)}
}
