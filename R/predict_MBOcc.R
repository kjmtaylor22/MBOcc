#' @title Predict from MBOcc object
#' @description Obtain 95% confidence intervals around predicted values
#' @author Kara J.M. Taylor (`k.taylor2@ufl.edu`)
#' @param MBOcc.obj.est Output list from `MBOcc` containing estimates and best-fit models
#' @param newdata Dataframe containing all covariates pertaining to MBOcc model
#' @param level Confidence interval 1-alpha
#' @param method One of 'default', 'delta', or 'boot', for different variance estimation methods
#' @export

predict.MBOcc.obj.est <- function(MBOcc.obj.est, newdata, level=0.95, each=F, method="default"){

  y <- lapply(MBOcc.obj.est$estimates, function(x){

    l <- length(unique(x$assigned.psi))

    form <- unique(x$call)

    if (length(form)>1){
      form <- paste(form, collapse=" ")
      form <- gsub(" ~1", "", form)
    }
    form <- terms(formula(form))
    mm <- model.matrix(form, newdata) #linear predictor matrix

    b <- x$mles #regression coefficient vector

    if (l>1){
      tmp <- mm
      i <- 1
      while (i < l){
        mm <- cbind(mm, tmp)
        i <- i+1
      }
    }

    if (each==T&l>1){
      yh <- sweep(mm, 2, b, FUN="*") #the fit with the new data
    } else {
      yh <- c(mm %*% b)
    }

    if (method=="boot"){
      if (form=="~1"){

        V <- MASS::ginv(x$hess) #variance-covariance matrix of mles
        if (each==T&l>1){
          yh <- c(yh)
          var.fit <- c(sweep(mm, 2, diag(V), FUN="*"))
        } else {
          var.fit <- c(mm %*% diag(V))
        }
        mu.y  <- 1/(1+exp(-(yh)))
        first.deriv <- exp(-mu.y)/(1+exp(-mu.y))^2

        var.p.x <- (first.deriv^2)*var.fit
        se.fit <- sqrt(var.p.x) #standard errors of the fitted mles

        alpha <- (1-level)/2 #Type I error rate for two-tailed distribution
        z <- qnorm(1-alpha) #Calculate appropriate z-value for given alpha

        lcl <- mu.y - se.fit*z #lower confidence interval around fitted model
        ucl <- mu.y + se.fit*z #upper confidence interval around fitted model

      } else {

        if (each==T&l>1){
          s <- split(1:length(b), rep(1:l, each=length(b)/l))
          tmp1 <- data.frame()
          tmp2 <- data.frame()
          for (i in 1:length(s)){
            tmp1 <- rbind(tmp1, mm[,s[[i]]])
            tmp2 <- rbind(tmp2, yh[,s[[i]]])
          }
          mm <- tmp1
          yh <- rowSums(tmp2)
        }

        n <- 200
        p  <- 1/(1+exp(-(yh)))

        boot <- matrix(NA, nrow(mm), n)
        for (i in 1:n){
          boot[,i] <- rbinom(length(p),1,p)
        }

        if (l>1&each==T){
          s <- split(1:nrow(boot), rep(1:l, each=nrow(boot)/l))
          try <- data.frame()
          for (i in 1:l){
            tmp <- apply(boot[s[[i]],], 2, function(x){
              glm.fit(mm[s[[i]],], x, family=binomial())$fitted.values
            })
            try <- rbind(try, tmp)
          }
        } else {
          try <- apply(boot, 2, function(x){
            glm.fit(mm, x, family=binomial())$fitted.values
          })
        }

        var.fit <- t(apply(try, 1, function(x){quantile(x, probs=c(0.025, 0.975))}))

        lcl <- var.fit[,1]
        ucl <- var.fit[,2]
        mu.y <- p

        alpha <- (1-level)/2 #Type I error rate for two-tailed distribution
        z <- qnorm(1-alpha) #Calculate appropriate z-value for given alpha
        se.fit <- (ucl - apply(try, 1, function(x){quantile(x, probs=c(0.5))}))/z
      }

    } else {

      V <- MASS::ginv(x$hess) #variance-covariance matrix of mles

      if (form=="~1"){
        if (each==T&l>1){
          yh <- c(yh)
          var.fit <- c(sweep(mm, 2, diag(V), FUN="*"))
        } else {
          var.fit <- c(mm %*% diag(V))
        }
      } else {
        var.fit <- (mm %*% V) * mm
        if (each==T&l>1){
          s <- split(1:length(b), rep(1:l, each=length(b)/l))
          tmp1 <- data.frame()
          tmp2 <- data.frame()
          for (i in 1:length(s)){
            tmp1 <- rbind(tmp1, var.fit[,s[[i]]])
            tmp2 <- rbind(tmp2, yh[,s[[i]]])
          }
          var.fit <- rowSums(tmp1)
          yh <- rowSums(tmp2)
        } else {var.fit <- rowSums(var.fit)} #variances for predicted means
      }

      alpha <- (1-level)/2 #Type I error rate for two-tailed distribution
      z <- qnorm(1-alpha) #Calculate appropriate z-value for given alpha

      if (method=="delta"){
        mu.y  <- 1/(1+exp(-(yh)))
        first.deriv <- exp(-mu.y)/(1+exp(-mu.y))^2
        var.fit <- (first.deriv^2)*var.fit
        se.fit <- sqrt(var.fit) #standard errors of the fitted mles

        lcl <- mu.y - se.fit*z #lower confidence interval around fitted model
        ucl <- mu.y + se.fit*z #upper confidence interval around fitted model

      } else {
        se.fit <- sqrt(var.fit) #standard errors of the fitted mles

        lcl <- yh - se.fit*z #lower confidence interval around fitted model
        ucl <- yh + se.fit*z #upper confidence interval around fitted model

        mu.y  <- 1/(1+exp(-(yh)))
        lcl  <- 1/(1+exp(-(lcl)))
        ucl  <- 1/(1+exp(-(ucl)))
      }
    }

    n <- nrow(newdata)
    newdata <- newdata[rep(1:n, l),]
    newdata <- cbind(newdata, Sites=rep(1:l, each=n))

    fit <- data.frame(formulae=paste(x$call, collapse=" "),
                      psi=paste(x$assigned.psi, collapse=" "),
                      newdata,
                      fit=mu.y,
                      se.fit=se.fit,
                      lcl=lcl,
                      ucl=ucl)

    out <- data.frame(formulae=paste(x$call, collapse=" "),
                      psi=paste(x$assigned.psi, collapse=" "),
                      betas=colnames(mm),
                      mles=b)

    return(fit)
  })

  class(y) <- "predMBOcc"
  return(y)
}
