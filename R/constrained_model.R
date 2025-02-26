#' @title Occupancy estimation for pathogen in four tissue types
#' @description
#' A short description...
#' @author Jose Miguel Ponciano Castellanos (josemi@ufl.edu)
#' @param formula of the form {~ Age * Species} to designate the covariate columns of interest
#' @param data the data, containing the covariates and the frequencies of state occupation
#' @param assign.psi a vector of the same length as the set of state positions (here CINT (4)) to indicate how the tissues should be grouped
#' @export

const.model <- function(formula, data, assign.psi){

  # formula: of the form {~ Age * Species} to designate the covariate columns of interest
  # data: the data, containing the covariates and the frequencies of state occupation
  # assign.psi: a vector of the same length as the set of state positions, here CINT (4)
  # to indicate how the tissues should be grouped

  mlogit<-function(a) {
    x=exp(a)/(1+sum(exp(a)))
    return(x)
  }

  S=length(assign.psi)        #  number of sub-levels of occupancy (4)
  SS=2^S               #  number of possible states, given sub-levels, including 0000
  SSm1=SS-1            #  number of possible states, given sub-levels excluding 0000

  #creating combination of truestates
  trustate=NULL
  for (i in 0:SSm1) trustate=rbind(trustate,as.integer(intToBits(i)[S:1]))
  states=CINT=apply(trustate,1,paste,collapse='') #CINT represent Cecum, Ileum, Nasal, Trachea


  ys.mat <- as.matrix(data[,grep("^[0-1]+$",names(data))])
  nstates <- ncol(ys.mat)
  nsamps  <- nrow(ys.mat)

  designmats <- list()
  nbetas.vec <- rep(0,S)

  for(i in 1:S){
    designmat <- model.matrix(formula[[i]], data=data)
    designmats[[i]] <- designmat
    nbetas.vec[[i]] <- ncol(designmat)
  }


  # Now tallying up the numper of parameters and stopping if >=nstates
  # Note: number of parameter sets == length of unique formulae BUT that still
  # leaves open the possibility of having two or more identical parameter sets
  # So the number of UNIQUE parameter sets is still given by:
  nparm.sets <- length(unique(assign.psi))
  betas.assign.df <- data.frame(betas=nbetas.vec,assign.psi=as.factor(assign.psi))
  levs.assign <- levels(betas.assign.df$assign.psi)
  distinct.betas <- rep(0,nparm.sets)
  for(i in 1:nparm.sets){
    subdat <- betas.assign.df[betas.assign.df$assign.psi == levs.assign[i],]
    distinct.betas[i] <- subdat[1,1]
  }
  nbetas <- sum(distinct.betas)

  # Total number of parameters
  if (nbetas >= nstates){
    warning("You have an overparameterized model. Reduce the number of covariates.")
    return(NULL)
  }

  # The length of this list of guesses is equal to the number of distinct sets of
  # parameters
  guesses.list <- list()
  for(i in 1:nparm.sets){
    guesses.list[[i]] <- runif(n=distinct.betas[i])
  }

  guesses<- unlist(guesses.list)

  #nll.calc <- function(guesses=guesses, SSm1=SSm1, ys.mat=ys.mat, trustate=trustate,
  #                     distinct.betas=distinct.betas, assign.psi=assign.psi,designmats=designmats){
  nll.calc <- function(guesses=guesses,out4optim=TRUE){


    nparm.sets <- length(distinct.betas)
    guesses.list2 <- list()
    i.start<-1
    for(i in 1:nparm.sets){
      i.end <- i.start+distinct.betas[i]-1
      guesses.list2[[i]] <- guesses[i.start:i.end]
      i.start <- i.end+1
    }

    # Sweeping vertically
    regs.output <- matrix(0,nrow=nsamps, ncol=S) # 4=because we have four tissues
    for(i in 1:S){
      designmat <- designmats[[i]]
      ith.betas  <- guesses.list2[[assign.psi[i]]]
      ith.output <- designmat%*%ith.betas
      regs.output[,i] <- ith.output
    }

    if (length(unique(assign.psi))!=length(assign.psi)){
      xs <- unlist(designmats[-which(duplicated(assign.psi))], recursive=F)
    } else {xs <- unlist(designmats, recursive=F)}
    xs.mat <- matrix(xs, nsamps, nbetas, byrow=F)


    # Sweeping horizontally but vertically (really)
    predpis.mat <- matrix(0,nrow=nsamps, ncol=SS)
    negllike.vec <- rep(0,nsamps)

    sweep <- function(x){

      psivec <- 1/(1+exp(-x))
      psivec[psivec>0.9999]  <- 0.99999
      psivec[psivec<0.00001] <- 0.00001

      # Create a matrix where I repeat that vector in each row, for 15 rows (16-1)
      # psi.mat <- matrix(rep(psivec,SSm1),nrow=SSm1,ncol=4,byrow=TRUE)
      psi.mat <- matrix(rep(psivec,SS),nrow=SS,ncol=S,byrow=TRUE)

      # Now create a `one minus psi` matrix where I repeat the vector of 1-psis 15 times
      # ompsi.mat <- matrix(rep(1-psivec,SSm1),nrow=SSm1,ncol=4,byrow=TRUE)
      ompsi.mat <- matrix(rep(1-psivec,SS),nrow=SS,ncol=S,byrow=TRUE)

      # Now elevate these matrices to the powers of the "trustate" and 1-"trustate"
      # pis.mat <- psi.mat^trustate[1:SSm1,]
      # onempis.mat <- ompsi.mat^(1-trustate[1:SSm1,])
      pis.mat <- psi.mat^trustate[1:SS,]
      onempis.mat <- ompsi.mat^(1-trustate[1:SS,])


      pis.prod <- pis.mat*onempis.mat
      pis.logit  <- apply(pis.prod,1,FUN=prod)

      ## And now we are ready to apply the inverse logit transform coded in the `mlogit' function:
      # pis.vec1 <- mlogit(pis.logit)
      ## which gives us a vector of length 15 which to which we then append the last value as the
      ## remainder from one

      #pis.vec <- c(pis.vec1,1-sum(pis.vec1))
      pis.vec1 <- pis.logit
      pis.vec <- pis.vec1/sum(pis.vec1)

      predpis.mat[i,] <- pis.vec
    }

    predpis.mat <- t(apply(regs.output, 1, sweep))

    for(i in 1:nsamps){

      ys <- ys.mat[i,]
      ith.N <- sum(ys)
      ith.negllike <- -dmultinom(x=ys, size=ith.N, prob=predpis.mat[i,], log=TRUE)

      negllike.vec[i] <- ith.negllike
    }

    tot.negllike <- sum(negllike.vec)
    #print(tot.negllike)
    if(out4optim==TRUE){
      return(tot.negllike)}else{

        outlist <- list(nll = tot.negllike, preds=predpis.mat, des=xs.mat, regs=regs.output)
        return(outlist)
      }
  }
  # Now we just optimize the negllike

  #optim.out <- optim(par=guesses, fn= nll.calc, method="BFGS",SSm1=SSm1, ys.mat=ys.mat, trustate=trustate,
  #                     distinct.betas=distinct.betas, assign.psi=assign.psi,designmats=designmats)


  optim.out <- optim(par=guesses, fn=nll.calc, method="BFGS", hessian=TRUE, out4optim=TRUE)

  call <- as.character(formula)
  assigned.psi <- assign.psi
  max.negll <- optim.out$value
  mles <- optim.out$par
  hess <- optim.out$hessian
  BIC <- 2*max.negll + length(mles)*log(nsamps)
  # add another function that returns the predictions
  get.mats <- nll.calc(guesses=mles, out4optim=FALSE)

  outbag <- list(call=call, assigned.psi=assigned.psi, max.negll=max.negll, mles=mles, BIC=BIC,
                 hess=hess, preds.mat=get.mats$preds, des.mat=get.mats$des, regs.mat=get.mats$regs)

  # Return a list with BIC, mles,max. negllike, predictions, etc...
  return(outbag)

}
