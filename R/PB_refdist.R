
### ###########################################################
###
### Parallel computing of reference distribution
###
### ###########################################################

PBrefdist <- function(largeModel, smallModel, nsim=200, seed=NULL, cl=NULL, details=0){
    UseMethod("PBrefdist")
}

PBrefdist.lm <- function(largeModel, smallModel, nsim=200, seed=NULL, cl=NULL, details=0){

  ##cat(".....PBrefdist.lm\n")
  t0 <- proc.time()
  .is.cluster <- !is.null(cl) && inherits(cl, "cluster")
    
  fun <- function(ll, ss, nsim, seed=seed){
    ##cat(".....fun\n")
    ## print(ss)
    ##     print(getCall(ss))
    ##     print(ss$data)
    
    simdata <- simulate(ss, nsim, seed=seed)
    ee  <- new.env()
    ee$simdata <- simdata

    cl.l <- getCall(ll)
    cl.s <- getCall(ss)

    ff.l <- update.formula(formula(ll),simdata[,ii]~.)
    ff.s <- update.formula(formula(ss),simdata[,ii]~.)
    environment(ff.l) <- environment(ff.s) <- ee 
    cl.l$formula <- ff.l
    cl.s$formula <- ff.s
    cl.l$data <- ss$data
    cl.s$data <- ss$data

    ref <- rep.int(NA,nsim)
    for (ii in 1:nsim){
      ref[ii] <- 2*(logLik(eval(cl.l))-logLik(eval(cl.s)))
    }
    ref    
  }
  
  if (!.is.cluster){
    ref <- fun(largeModel, smallModel, nsim, seed=seed)
  } else {
    nsim2 <- round(nsim/length(cl))
    if (details>=1)
      cat(sprintf("* Using %i clusters and %i samples per cluster\n", length(cl), nsim2))
    clusterExport(cl, ls(envir=.GlobalEnv), envir = .GlobalEnv)
    #clusterSetupSPRNG(cl)
    clusterSetRNGStream(cl)
    ref <- unlist(clusterCall(cl, fun, largeModel, smallModel, nsim2))
  }

  ref <- ref[ref>0]
  ctime <- (proc.time()-t0)[3]
  attr(ref,"ctime") <- ctime
  if (details>0)
    cat(sprintf("Reference distribution with %i samples; computing time: %5.2f secs. \n",
                length(ref), ctime))
  
  ref
}

PBrefdist.mer <- function(largeModel, smallModel, nsim=200, seed=NULL, cl=NULL, details=0){

    t0 <- proc.time()

    ##if(is.null(smallModel@call$REML) || smallModel@call$REML)
    if (getME(smallModel, "is_REML"))
      smallModel <- update(smallModel,REML=FALSE)

    ##if(is.null(largeModel@call$REML) || largeModel@call$REML)
    if (getME(largeModel, "is_REML"))
      largeModel <- update(largeModel,REML=FALSE)

    .is.cluster <- !is.null(cl) && inherits(cl, "cluster")
    
    fun <- function(lg, sm, nsim=200, seed=NULL){
      simdata <- as.matrix(simulate(sm, nsim=nsim, seed=seed))
      ref     <- rep(NA, nsim)
      for (kk in seq_len(nsim)){
        yyy   <- simdata[,kk]
##         small <- suppressWarnings(refit(ss, newresp=yyy))
##         large <- suppressWarnings(refit(ll, newresp=yyy))
        sm2 <- refit(sm, newresp=yyy)
        lg2 <- refit(lg, newresp=yyy)
        ttt   <- as.numeric(2*(logLik(lg2, REML=FALSE) -
                               logLik(sm2, REML=FALSE)))
        ref[kk]  <- ttt
      }
      ref
    }
    
    if (!.is.cluster){
      ref <- fun(largeModel, smallModel, nsim, seed=seed)
    } else {      
      nsim2 <- round(nsim/length(cl))
      if (details>=1)
        cat(sprintf("* Using %i clusters and %i samples per cluster\n", length(cl), nsim2))
      clusterEvalQ(cl, library(lme4))
      clusterSetRNGStream(cl)
      ref <- unlist(clusterCall(cl, fun, largeModel, smallModel, nsim2))
    }

    neg <- sum(ref<0)
##     cat(sprintf("neg=%f\n", neg/nsim))
    
    ref   <- ref[ref>0]
    ctime <- (proc.time()-t0)[3]
    attr(ref,"ctime") <- ctime
    if (details>0)
      cat(sprintf("Reference distribution with %i samples; computing time: %5.2f secs. \n",
                  length(ref), ctime))
    
    ref
}

