#' @title Predict from MBOcc object
#' @description Obtain 95% confidence intervals around predicted values
#' @author Kara J.M. Taylor (`k.taylor2@ufl.edu`)
#' @param MBOcc.obj.est Output list from `MBOcc` containing estimates and best-fit models
#' @param newdata Dataframe containing all covariates pertaining to MBOcc model
#' @param level Confidence interval 1-alpha
#' @export

predict.MBOcc.obj.est <- function(MBOcc.obj.est, newdata, level=0.95){

  y <- lapply(MBOcc.obj.est$estimates, function(x){

    form <- unique(x$call)
    if (length(form)>1){
      form <- paste(form, collapse=" ")
      form <- gsub(" ~1", "", form)
    }
    form <- terms(formula(form))
    mm <- model.matrix(form, newdata) #linear predictor matrix
    b <- x$mles #regression coefficient vector
    yh <- c(mm %*% b) #the fit with the new data

    V <- matlib::inv(x$hess) #variance-covariance matrix of mles
    var.fit <- rowSums((mm %*% V) * mm) #variances for predicted means
    se.fit <- sqrt(var.fit) #standard errors of the fitted mles

    alpha <- (1-level)/2 #Type I error rate for two-tailed distribution
    z <- qnorm(1-alpha) #Calculate appropriate z-value for given alpha

    lcl <- yh - se.fit*z #lower confidence interval around fitted model
    ucl <- yh + se.fit*z #upper confidence interval around fitted model

    fit <- data.frame(formulae=paste(x$call, collapse=" "),
                      psi=paste(x$assigned.psi, collapse=" "),
                      newdata,
                      fit=yh,
                      se.fit=se.fit,
                      lcl=lcl,
                      ucl=ucl)

    out <- data.frame(formulae=paste(x$call, collapse=" "),
                      psi=paste(x$assigned.psi, collapse=" "),
                      betas=colnames(mm),
                      mles=b)

    return(fit)
  })

  return(y)
}
