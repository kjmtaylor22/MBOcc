#' @title Estimate time to model completion
#' @param n.sp The number of species
#' @param n.form The number of formulae
#' @param n.psi The number of unique psi assignments
#' @export

MBruntime <- function(n.sp, n.form, n.psi){

  sp <- log10(c(2,16,60,361,2064))
  newdata <- expand.grid(psi=c(1:5),
                       form=c(2,4,6),
                       sp=sp)

  train <- data.frame(newdata,
                        runtime=log10(
                         c(0.805,2.78,6.063,9.764,14.95,9.501,17.027,29.555,39.689,
                           43.487,17.749,27.518,47.987,67.399,72.414,5.939,20.242,
                           49.263,76.991,121.885,87.128,142.669,238.165,331.573,
                           380.358,137.593,219.523,397.072,565.97,611.179,22.047,
                           74.614,178.325,278.9,437.789,247.234,417.854,728.54,
                           1041.818,1251.531,390.871,700.187,1228.689,1810.935,
                           1984.01,98.144,322.868,692.772,1079.517,1670.174,942.015,
                           1505.443,2662.768,3821.902,4500.767,1414.023,2382.089,
                           4320.619,6293.292,6870.837,321.706,1134.647,2711.463,
                           4157.906,6675.549,4577.276,7227.623,12797.01,18742.85,
                           20641.35,6907.341,11416.35,20945.35,30983.38,33640.44)/60))

  #par(xpd=T)
  #base::plot(train$psi+train$form, train$runtime, col=rep(1:5, each=15), pch=20)
  #legend("top", inset=c(0.5,-0.4), horiz=T, legend=as.character(10^sp), col=1:5, pch=20, title="n Species")

  ## get y-intercepts and slopes from each curve in above plot
  intercepts <- c()
  slopes <- c()
  asym <- c()
  for (i in sp){
    #m <- lm(runtime ~ log(I(psi+form)), data=subset(train, subset=sp==i))
    m <- nls(runtime ~ SSasymp(I(psi+form), Asym, R0, lrc), data=subset(train, subset=sp==i))
    intercepts <- c(intercepts, coef(m)[2])
    slopes <- c(slopes, coef(m)[3])
    asym <- c(asym, coef(m)[1])
    #curve(coef(m)[2] * log(x) + coef(m)[1], add=T, col=which(sp==i))
    #curve(coef(m)[1]-(coef(m)[1]-coef(m)[2]) * exp(-exp(coef(m)[3])*x), add=T, col=which(sp==i))
  }

  ## see how the intercepts change with the number of species in the dataset
  #plot(sp, intercepts, pch=0, cex=1.5)


  ## define candidate models
  #pred1 <- nls(intercepts ~ a/(1 + exp(-b * (sp))), start=list(a=-1,b=.005))
  #pred2 <- nls(intercepts ~ a * (sp + b) ^ (1/3) - c, start = list(a=1, b=0, c=0))
  pred3 <- nls(intercepts ~ c+(d-c)/(1+exp(b*(sp-1.56))), start = list(d=-2,c=-5,b=1))
  pred4 <- lm(asym ~ sp)

  ## add the predictions to the graph
  #points(sp, predict(pred1), col="red", pch=2)
  #points(sp, predict(pred2), col="blue", pch=6)
  #points(sp, predict(pred3), col="magenta", pch=4)
  #curve(coef(pred1)[1]/(1 + exp(-coef(pred1)[2] * (x))) , add=T, col="red", lwd=2)
  #curve(coef(pred2)[1] * (x + coef(pred2)[2]) ^ (1/3) - coef(pred2)[3], add=T, col="blue", lwd=2)
  #curve(coef(pred3)[2]+(coef(pred3)[1]-coef(pred3)[2])/(1+exp(coef(pred3)[3]*(x-1.56))), add=T, col="magenta", lwd=2)
  #legend("top", inset=c(0.5,-0.4), horiz=T, legend=c("Observed","Logistic","Cube root"),
  #       col=c("black","red","blue"), pch=c(0,NA,NA), lty=c(0,1,1), lwd=c(0,2,2), title="Observed vs Predicted")

  ## compare candidate models
  #anova(pred1, pred2) # not significantly different, both pretty good

  ## Calculate weights for the models
  #delta <- BIC(pred2) - BIC(pred1)
  #weight1 <- exp(-0.5 * delta) / (1 + exp(-0.5 * delta))
  #weight2 <- 1 - weight1

  ## use weights to calculate a weighted y-intercept and y-shift for the user n.sp parameter
  #yint <- predict(pred1, newdata=data.frame(sp=0))*weight1+predict(pred2, newdata=data.frame(sp=0))*weight2
  #shift <- predict(pred1, newdata=data.frame(sp=log10(n.sp)))*weight1+
  #  predict(pred2, newdata=data.frame(sp=log10(n.sp)))*weight2-yint
  r0 <- predict(pred3, newdata=data.frame(sp=0))
  lrc <- mean(slopes)
  shift <- predict(pred3, newdata=data.frame(sp=log10(n.sp))) - r0

  ## generate a baseline for n.sp=1 (log(n.sp)=0) from which to shift
  baseline <- expand.grid(psi=c(1:5),form=c(2,4,6),sp=0)
  #baseline <- dplyr::mutate(baseline, runtime=yint+mean(slopes)*log(psi+form))
  baseline <- dplyr::mutate(baseline,
                            runtime=coef(pred4)[1]-(coef(pred4)[1]-r0) * exp(-exp(lrc)*(psi+form)))

  ## generate baseline model from baseline data
  #try <- lm(runtime ~ log(I(psi+form)), data=baseline)

  ## this is the user input parameters
  #testset <- data.frame(form=1:n.form, psi=n.psi)

  ## predict runtime based on baseline model and n.sp-dependent y-shift
  #runtime <- predict(try, newdata = testset) + shift
  a <- coef(pred4)[1]+coef(pred4)[2]*log10(n.sp)
  runtime <- a-(a-r0) * exp(-exp(lrc)*(n.psi+n.form))

  ## How does this fit on another visualization of the runtime data?
  #plot(train$psi+train$form, train$runtime, col=rep(1:5, each=15), pch=20, xlim=c(0,15), ylim=c(-3,3))
  #points(baseline$form+baseline$psi, baseline$runtime, pch=16, col="grey40")
  #points(n.form+n.psi, runtime, pch=8, col="grey40")

  ## Convert estimated runtime to minutes
  estimate <- 10^runtime
  t.un <- "minutes"

  ## Convert from minutes to hours or seconds if necessary
  if(estimate>60){
    estimate=estimate/60; t.un="hours"
  }
  if(estimate<1){
    estimate=estimate*60; t.un="seconds"
  }

  ## output console message to the user
  message(paste("Estimated runtime: ", round(estimate, 1), t.un))
}
