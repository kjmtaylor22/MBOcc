#' @title Plot MBOcc object
#' @description Plotting function for MBOcc object
#' @author Kara J.M. Taylor (`k.taylor2@ufl.edu`)
#' @param MBOcc MBOcc object
#' @param plot Which parameter estimates to plot, one of `MLE`, `Site`, or `State`
#' @import ggplot2
#' @import dplyr
#' @import ggthemes
#' @export

plot <- function(MBOcc, plot=c("MLE","Site","State"), return=F){

  library(ggplot2)

  stopifnot(class(MBOcc)=="MBOcc", plot%in%c("MLE","Site","State"))

  colsSP <- group.colors(unique(as.character(MBOcc$P.site$L1)), show=F)

  groups <- MBOcc$groups

  if (return==T){out <- list()}

  if (plot=="MLE"){

    if (return==F){
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }

    mle0 <- MBOcc$MLE.SE

    plots <- unique(mle0[,c("model.psi","site.psi","site.call")])
    for (i in 1:nrow(plots)){
      sub <- droplevels(subset(mle0, subset=model.psi==plots[i,1]&site.psi==plots[i,2]&site.call==plots[i,3]))

      g <- ggplot(sub) +
        labs(title=paste(paste0("Assigned psi: ", plots[i,1]),
                         paste0("Site psi: ", plots[i,2]),
                         paste0("Call: ", plots[i,3]), sep="\n")) +
        facet_wrap(vars(beta)) +
        geom_hline(yintercept=0.5, color="grey60", linetype="dashed") +
        geom_segment(aes(x=L1, xend=L1, y=-Inf, yend=exp(mle)/(1+exp(mle))), color="grey80", linetype="dotted") +
        geom_point(aes(L1, exp(mle)/(1+exp(mle)), color=L1), show.legend = F) +
        geom_errorbar(aes(x=L1, ymin=exp(lcl)/(1+exp(lcl)), ymax=exp(ucl)/(1+exp(ucl)), color=L1), show.legend = F) +
        geom_text(aes(x=L1, y=exp(ucl)/(1+exp(ucl)), label=sig), color="black", size=5, nudge_y=0.01) +
        ylab(expression("MLE of " ~beta)) + xlab("") +
        scale_color_manual(values=colsSP) + ylim(0,1) + coord_cartesian(ylim=c(-0.01,1.01), expand = F) +
        ggthemes::theme_few() +
        theme(axis.text.x = element_text(angle=75, vjust=1, hjust=1),
              title=element_text(size=8))

      if (return==T){out[[paste(i, plots[i,3])]] <- g} else {print(g)}
    }
  }

  if (plot=="Site"){

    library(dplyr)

    if (return==F){
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }

    site0 <- MBOcc$P.site

    plots <- unique(site0[,c("variable","formulae","psi")])

    forms <- strsplit(plots$formulae, split="~")
    forms <- lapply(forms, function(x){gsub(",| ","",x)[-1]})
    forms <- unlist(forms[-which(duplicated(plots[,c(2,3)]))])

    terms <- strsplit(plots$psi, split=".", fixed = T)
    terms <- unlist(terms[-which(duplicated(plots[,c(2,3)]))])

    plots <- cbind(plots, call=forms, site.psi=terms) %>%
      group_by(psi, call, site.psi) %>%
      mutate(site=paste(variable, collapse="/")) %>%
      as.data.frame()

    plots <- unique(plots[,-1])


    for (i in 1:nrow(plots)){
      sub <- droplevels(subset(site0, subset=psi==plots[i,2]&formulae==plots[i,1]))

      for (j in groups){
        pull <- grep(j, plots$call[i])

        if (any(!is.na(pull))){
          sub.pull <- sub
          names(sub.pull)[which(names(sub.pull)==j)] <- "covariate"

          p1 <- as.numeric(plots[i,4])
          p2 <- plots[i,5]

          g1 <- ggplot(sub.pull) +
            labs(title=paste(paste("Assigned psi:", plots[i,2]), paste("Call:", plots[i,3]), sep="\n"),
                 subtitle=substitute(paste(psi[p1], ": ", p2))) +
            facet_wrap(vars(L1)) +
            geom_hline(yintercept=0.5, color="grey80", linetype="dotted") +
            geom_point(aes(covariate, exp(value)/(1+exp(value)), color=L1), shape=20, size=2,
                       position=position_jitter(width=0.2), show.legend = F) +
            ylab(expression("Probability of occupancy " ~psi)) +
            xlab(j) +
            scale_color_manual(values=colsSP) +
            ggthemes::theme_few() +
            theme(title=element_text(size=8))

          if (return==T){out[[paste(i, j)]] <- g1} else {print(g1)}
        }
      }
    }

  }

  if (plot=="State"){

    if (return==F){
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }

    state0 <- MBOcc$P.state

    foo <- unique(state0[,c("formulae","psi")])

    forms <- strsplit(foo$formulae, split="~")
    forms <- lapply(forms, function(x){gsub(",| ","",x)[-1]})

    terms <- strsplit(foo$psi, split=".", fixed = T)

    plots <- data.frame()
    for (x in 1:length(terms)){
      tmp <- suppressWarnings(cbind(foo[x,], call=forms[[x]], site.psi=terms[[x]]))
      tmp <- unique(tmp)
      plots <- rbind(plots, tmp)
    }

    plots <- plots[-which(duplicated(plots[,1:3])),]

    for (i in 1:nrow(plots)){
      sub <- droplevels(subset(state0, subset=psi==plots[i,2]&formulae==plots[i,1]))

      for (j in groups){
        pull <- grep(j, plots$call[i])

        if (any(!is.na(pull))){
          sub.pull <- sub
          names(sub.pull)[which(names(sub.pull)==j)] <- "covariate"

          if (class(sub.pull$covariate)%in%c("factor","character","logical")){
            geom <- expression(stat_summary(aes(as.factor(covariate), value , fill=variable),
                                 fun="mean", geom="col", position="stack"))
          } else {
            geom <- expression(stat_summary(aes(covariate, value, fill=variable),
                                 color="white", linewidth=0.03,
                                 fun="mean", geom="area", position="stack"))
          }

          p2 <- plots[i,3]

          g2 <- ggplot(sub.pull) +
            labs(title=paste(paste0("Assigned psi: ", plots[i,2])),
                 subtitle=substitute(paste("Call: ", psi, " ~ ", p2))) +
            facet_wrap(vars(L1)) +
            geom_hline(yintercept=0, color="grey80", linetype="dotted") +
            eval(geom) +
            ylab(expression("State probability " ~pi)) +
            xlab(j) +
            scale_fill_viridis_d(option="magma", name="State") +
            ggthemes::theme_few() +
            theme(title=element_text(size=8))

          if (return==T){out[[paste(i, j)]] <- g2} else {print(g2)}
        }
      }
    }

    if (any(plots$call==1)){

      library(dplyr)

      pull <- plots[plots$call==1,]
      sub <- subset(state0, subset=formulae%in%pull$formulae&psi%in%pull$psi)

      for (i in unique(pull$psi)){
        sub.pull <- subset(sub, subset=psi==i)

        g2 <- ggplot(sub.pull) +
          labs(title=paste(paste0("Assigned psi: ", i)),
               subtitle=substitute(paste("Call: ", psi, " ~ 1 (intercept-only)"))) +
          geom_hline(yintercept=0, color="grey80", linetype="dotted") +
          stat_summary(aes(L1, value, fill=variable), fun="mean", geom="col", position="stack") +
          ylab(expression("State probability " ~pi)) +
          xlab("Species") +
          scale_fill_viridis_d(option="magma", name="State") +
          ggthemes::theme_few() +
          theme(axis.text.x = element_text(angle=75, vjust=1, hjust=1),
                title=element_text(size=8))

        if (return==T){out[[i]] <- g2} else {return(g2)}
      }
    }
  }
  if (return==T){return(out)}
}
