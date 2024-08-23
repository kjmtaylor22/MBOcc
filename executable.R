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
groups <- c("Experiment","Species","Age_weeks","Flock")

test <- format(comm, meta, tax, ids, groups, 1000, "BodySite")


formula <- list(c(~1,~1,~1,~1),
                #c(~1+Age_weeks,~1+Age_weeks,~1+Age_weeks,~1+Age_weeks),
                #c(~1+Flock,~1+Flock,~1+Flock,~1+Flock))#,
                c(~1+Age_weeks*Flock,~1+Age_weeks*Flock,~1+Age_weeks*Flock,~1+Age_weeks*Flock),
                c(~1+Age_weeks*Species,~1+Age_weeks*Species,~1+Age_weeks*Species,~1+Age_weeks*Species))

formulas <- list(formula)#, formula, formula, formula, formula)
assigns <- list(c(1,1,1,1))#, c(1,1,2,2), c(1,2,3,3), c(1,1,2,3), c(1,2,3,4))

new <- MBOcc(test, formulas, assigns)

bag <- unroll(test, new, groups, as.character(unique(meta$BodySite)))

save(test, formulas, assigns, new, bag, file="data/MBOcc_output.RD")

load("data/MBOcc_output.RD")


plot(bag, plot="MLE")
plot(bag, plot="Site")
plot(bag, plot="State")

newdata <- data.frame(Age_weeks=rep(1:100, 6), Flock=rep(c("F1","F2","SPF-T","RN","RS","SPF-C"), each=100), Species=rep(c("Turkey","Chicken"), each=300))
try <- predict(new, newdata)
preds <- reshape2::melt(try, id.vars=c("Flock","Species","Age_weeks","formulae","psi")) %>% reshape2::dcast(formulae+psi+L1+Flock+Species+Age_weeks ~ variable)
colsSP <- group.colors(unique(as.character(bag$P.site$L1)))

g <- ggplot(preds) +
  facet_wrap(vars(formulae, L1)) +
  geom_ribbon(aes(x=Age_weeks, ymin=exp(lcl)/(1+exp(lcl)), ymax=exp(ucl)/(1+exp(ucl)), fill=Flock), alpha=0.3) +
  geom_smooth(aes(x=Age_weeks, y=exp(fit)/(1+exp(fit)), color=Flock), linewidth=1)

