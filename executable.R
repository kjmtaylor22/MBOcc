## Load the data

load("data/basefiles_26-03-24.RD")

library(dplyr)
library(ggplot2)

## Do some data cleanup

meta$Flock <- factor(meta$Flock, labels=c("","F1","F1","F2","F2","R","RN","RS","SPF-C","SPF-T","UMN"))

meta <- subset(meta, subset=Experiment!="TK-107")
meta <- mutate(meta, Rearing=NA)
meta$Rearing[grep("-85", meta$Experiment)] <- "Commercial"
meta$Rearing[-grep("-85", meta$Experiment)] <- "Research"

meta <- subset(meta, subset=BodySite%in%c("CECUM","ILEUM","NASAL","TRACHEA")&Species!="Control")

nas.expand <- meta[grep("-", meta$Bird_ID),]
nas.expand$Bird_ID <- gsub("--", "-", nas.expand$Bird_ID)
nas.expand <- tidyr::separate(nas.expand, "Bird_ID", c("from","to"),"-", convert=T)

out.expand <- data.frame()
for (j in 1:nrow(nas.expand)){
  k <- c(nas.expand$from[j]:nas.expand$to[j])

  out1 <- data.frame()
  for (l in 1:length(k)){out1 <- rbind(out1, nas.expand[j,])}

  out1$from <- k
  out1 <- out1[,-which(names(out1)=="to")]

  out.expand <- rbind(out.expand, out1)
}
names(out.expand)[names(out.expand)=="from"] <- "Bird_ID"

comm.expand <- data.frame()
for (i in out.expand$SampleID){
  comm.expand <- rbind(comm.expand, comm[i,])
}

comm <- comm[-match(nas.expand$SampleID, row.names(comm)),]
comm <- rbind(comm, comm.expand)

out.expand$SampleID <- row.names(comm.expand)
row.names(out.expand) <- out.expand$SampleID

meta <- meta[-grep("-", meta$Bird_ID),]
meta <- rbind(meta, out.expand)

metaR <- subset(meta, subset=Flock=="R")
metaR <- subset(metaR, subset=Age!="01W")
metaRN <- subset(metaR, subset=Bird_ID%in%c(41:50))
metaRN$Flock <- "RN"
metaRS <- subset(metaR, subset=Bird_ID%in%c(51:60))
metaRS$Flock <- "RS"

meta <- subset(meta, subset=Flock!="R")
meta <- rbind(meta, metaRN, metaRS)

row.names(meta) <- meta$SampleID
row.names(tax) <- tax$tag
comm <- comm[row.names(meta),]


## Actually run the models now

library(MBOcc)

ids <- c("Bird_ID")
groups <- c("Experiment","Species","Age_weeks","Flock","Rearing")

test <- format(comm, meta, tax, ids, groups, 1000, "BodySite")


formula <- list(c(~1,~1,~1,~1),
                #c(~1+Age_weeks,~1+Age_weeks,~1+Age_weeks,~1+Age_weeks),
                #c(~1+Flock,~1+Flock,~1+Flock,~1+Flock))#,
                c(~1+Age_weeks*Flock,~1+Age_weeks*Flock,~1+Age_weeks*Flock,~1+Age_weeks*Flock),
                c(~1+Age_weeks*Species,~1+Age_weeks*Species,~1+Age_weeks*Species,~1+Age_weeks*Species),
                #c(~1+Age_weeks*Experiment,~1+Age_weeks*Experiment,~1+Age_weeks*Experiment,~1+Age_weeks*Experiment),
                c(~1+Age_weeks*Rearing,~1+Age_weeks*Rearing,~1+Age_weeks*Rearing,~1+Age_weeks*Rearing))

formulas <- list(formula, formula, formula, formula, formula)
assigns <- list(c(1,1,1,1), c(1,1,2,2), c(1,2,3,3), c(1,1,2,3), c(1,2,3,4))

new <- MBOcc(test, formulas, assigns)

bag <- unroll(test, new, groups, as.character(unique(meta$BodySite)))

save(test, formulas, assigns, new, bag, file="output/MBOcc_output_full.RD")

load("output/MBOcc_output_full.RD")


mleplots <- plot(bag, plot="MLE", return=T)
siteplots <- plot(bag, plot="Site", return=T)
stateplots <- plot(bag, plot="State", return=T)

newdata <- data.frame(Age_weeks=rep(1:100, 6),
                      Experiment=rep(c("TK-85","TK-85","TK-93","CK-85","CK-85","CK-106"), each=100),
                      Species=rep(c("Turkey","Chicken"), each=300),
                      Rearing=rep(rep(c("Commercial","Commercial","Research"), each=100), 2),
                      Flock=rep(c("F1","F2","SPF-T","RN","RS","SPF-C"), each=100))
try <- predict(new, newdata, each = T)
try <- lapply(try, function(x){if (!"Sites"%in%names(x)){cbind(x,Sites=rep(1,nrow(x)))}else{return(x)}})
preds <- reshape2::melt(try, id.vars=c("Experiment","Species","Flock","Rearing","Age_weeks","formulae","psi","Sites")) %>%
  unique() %>%
  reshape2::dcast(formulae+psi+L1+Experiment+Species+Flock+Rearing+Age_weeks+Sites ~ variable)
colsSP <- group.colors(unique(as.character(bag$P.site$L1)))

g <- ggplot(preds) +
  facet_wrap(vars(psi, Sites, L1)) +
  geom_ribbon(aes(x=Age_weeks, ymin=exp(lcl)/(1+exp(lcl)), ymax=exp(ucl)/(1+exp(ucl)), fill=Experiment), alpha=0.3) +
  geom_smooth(aes(x=Age_weeks, y=exp(fit)/(1+exp(fit)), color=Experiment), linewidth=1) +
  paletteer::scale_color_paletteer_d("yarrr::appletv") + #see https://r-graph-gallery.com/color-palette-finder for colors
  paletteer::scale_fill_paletteer_d("yarrr::appletv") +
  ggthemes::theme_few()

save(mleplots, siteplots, stateplots, g, preds, file="data/MBOcc_plots_full.RD")

g1 <- egg::ggarrange(mleplots[[11]]+theme(axis.text.x=element_blank()),
                     mleplots[[12]]+theme(axis.text.x=element_blank()),
                     mleplots[[13]],mleplots[[14]])

g2 <- stateplots[[1]]+
  facet_grid(L1 ~ Flock, scales="free_x")  +
  theme(strip.text.y = element_text(size=7))
g3 <- egg::ggarrange(
  stateplots[[3]]+
    facet_grid(L1 ~ Rearing, scales="free_x")  +
    guides(fill="none") +
    theme(strip.text.y = element_text(size=7)),
  stateplots[[6]]+
    facet_grid(L1 ~ Species, scales="free_x")  +
    guides(fill="none") +
    theme(strip.text.y = element_text(size=7)),
  stateplots[[7]]+
    facet_grid(L1 ~ Rearing, scales="free_x")  +
    theme(strip.text.y = element_text(size=7)),
  nrow=1)

g4 <- ggplot(subset(preds, subset=psi=="1 1 1 1")) +
    facet_wrap(vars(L1)) +
    geom_ribbon(aes(x=Age_weeks, ymin=exp(lcl)/(1+exp(lcl)), ymax=exp(ucl)/(1+exp(ucl)), fill=Flock), alpha=0.3) +
    geom_line(aes(x=Age_weeks, y=exp(fit)/(1+exp(fit)), color=Flock), linewidth=1) +
    paletteer::scale_color_paletteer_d("yarrr::appletv") + #see https://r-graph-gallery.com/color-palette-finder for colors
    paletteer::scale_fill_paletteer_d("yarrr::appletv") +
    ggthemes::theme_few()
g5 <- ggplot(subset(preds[grep("Species",preds$formulae),], subset=psi=="1 2 3 3")) +
    facet_grid(Sites ~ L1) +
    geom_ribbon(aes(x=Age_weeks, ymin=exp(lcl)/(1+exp(lcl)), ymax=exp(ucl)/(1+exp(ucl)), fill=Species), alpha=0.3) +
    geom_smooth(aes(x=Age_weeks, y=exp(fit)/(1+exp(fit)), color=Species), linewidth=1) +
    paletteer::scale_color_paletteer_d("yarrr::appletv") + #see https://r-graph-gallery.com/color-palette-finder for colors
    paletteer::scale_fill_paletteer_d("yarrr::appletv") +
    ggthemes::theme_few()
g6 <- ggplot(subset(preds[grep("Rearing",preds$formulae),], subset=psi=="1 2 3 3")) +
  facet_grid(Sites ~ L1) +
  geom_ribbon(aes(x=Age_weeks, ymin=exp(lcl)/(1+exp(lcl)), ymax=exp(ucl)/(1+exp(ucl)), fill=Rearing), alpha=0.3) +
  geom_smooth(aes(x=Age_weeks, y=exp(fit)/(1+exp(fit)), color=Rearing), linewidth=1) +
  paletteer::scale_color_paletteer_d("yarrr::appletv") + #see https://r-graph-gallery.com/color-palette-finder for colors
  paletteer::scale_fill_paletteer_d("yarrr::appletv") +
  ggthemes::theme_few()
g7 <- ggplot(subset(preds, subset=psi=="1 1 2 3")) +
  facet_grid(Sites ~ L1) +
  geom_ribbon(aes(x=Age_weeks, ymin=exp(lcl)/(1+exp(lcl)), ymax=exp(ucl)/(1+exp(ucl)), fill=Rearing), alpha=0.3) +
  geom_smooth(aes(x=Age_weeks, y=exp(fit)/(1+exp(fit)), color=Rearing), linewidth=1) +
  paletteer::scale_color_paletteer_d("yarrr::appletv") + #see https://r-graph-gallery.com/color-palette-finder for colors
  paletteer::scale_fill_paletteer_d("yarrr::appletv") +
  ggthemes::theme_few()
