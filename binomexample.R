

# Simulating data 'on the go':


ntrials <- 1; # Bernoulli trials
nreps <- 200;
x <- runif(n=200, min=-3, max=3)
hist(x)

# Setting P(success) as a function of a covariate
beta0<- 1.5
beta1 <- 2.85;
real.p <- 1/(1+exp(-(beta0+beta1*x)))
plot(x, real.p, pch=16)

# Simulating data
data.sim <- rbinom(n=nreps, size=ntrials, prob=real.p)

# raw data
my.data <- cbind(x,data.sim)
colnames(my.data) <- c("covariate", "Successes")

# Negative log-likelihood function:
negll.logist <- function(guess, obs.vec, xvals){
	
	beta0 <- guess[1]
	beta1 <- guess[2]
	p.x  <-  1/(1+exp(-(beta0+beta1*x)))
	
	llike <- dbinom(x=obs.vec, size=1, prob=p.x, log=TRUE)
	negll <- -sum(llike)
	return(negll)
}

# Compute the neg log likelihood for a vectore of guesses of the parameter values
# This is just to see if the negloglike function "negll.logist" works properly
my.guess <- c(0.44,0.6)
negll.logist(guess=my.guess, obs.vec=my.data[,2], xvals=my.data[,1])

# Now the optimization itself
# Work with the neg-loglikelihood and then, you don't need to multiply the hessian by -1

ml.estim <- optim(par=my.guess, fn=negll.logist, method="Nelder-Mead", obs.vec=my.data[,2], xvals= my.data[,1], hessian=TRUE)

mles <- ml.estim$par
my.hess <- ml.estim$hessian

library("MASS") # for 'ginv'

Fish.Inv <- ginv(my.hess) #solve(my.hess)

zalphahalf <- qnorm(p=0.975, mean=0,sd=1) 
st.errs <- zalphahalf*sqrt(diag(Fish.Inv))

low.cis <- mles - st.errs
hi.cis  <- mles + st.errs

CIs.mat <- cbind(low.cis, mles, hi.cis)
colnames(CIs.mat)  <- c("2.5%", "MLE", "97.5%")
row.names(CIs.mat) <- c("Beta0", "Beta1")
print(CIs.mat)


beta0.hat <- mles[1]
beta1.hat <- mles[2]

nice.xs <- sort(x)
p.hat <- 1/(1+exp(-(beta0.hat +beta1.hat*nice.xs)))


mu.y  <- 1/(1+exp(-(beta0.hat +beta1.hat*nice.xs) )) 
var.y <- Fish.Inv[1,1] + (nice.xs^2)*Fish.Inv[2,2] + 2*nice.xs*Fish.Inv[1,2]

first.deriv <- exp(-mu.y)/(1+exp(-mu.y))^2

var.p.x <- (first.deriv^2)*var.y

L.CI <- p.hat - zalphahalf*sqrt(var.p.x)
U.CI <- p.hat + zalphahalf*sqrt(var.p.x)


par(oma=c(1,1,1,1), mar=c(4,4,2,1))
plot(nice.xs, p.hat, type="l", col="red", lwd=2, ylim=c(0,1.75), bty="l", 
xlab="x", ylab="p(x)", cex.lab=1.75)
points(nice.xs, L.CI, type="l", col="blue")
points(nice.xs, U.CI, type="l", col="blue")
abline(h=1, lty=2, col="black")
abline(h=0, lty=2, col="black")

 