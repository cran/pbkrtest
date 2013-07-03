
### ###########################################################
###
### Parallel computing of reference distribution
###
### ###########################################################

PBrefdist <- function(largeModel, smallModel, nsim=1000, seed=NULL, cl=NULL, details=0){
    UseMethod("PBrefdist")
}

PBrefdist.lm <- function(largeModel, smallModel, nsim=1000, seed=NULL, cl=NULL, details=0){

  ##cat(".....PBrefdist.lm\n")
  t0 <- proc.time()
  .is.cluster <- !is.null(cl) && inherits(cl, "cluster")
    
  fun <- function(ll, ss, nsim, seed=seed){    
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


.getref <- function(lg, sm, nsim=20, seed=NULL){
  simdata <- as.matrix(simulate(sm, nsim=nsim, seed=seed))
  ref     <- rep(NA, nsim)
  for (kk in seq_len(nsim)){
    yyy   <- simdata[,kk]
    sm2 <- refit(sm, newresp=yyy)
    lg2 <- refit(lg, newresp=yyy)
    ttt   <- 2*(logLik(lg2, REML=FALSE) - logLik(sm2, REML=FALSE))
    ref[kk]  <- ttt
  }
  as.numeric(ref)
}

PBrefdist.lmerMod <- PBrefdist.mer <- function(largeModel, smallModel, nsim=1000, seed=NULL, cl=NULL,details=0){

  t0 <- proc.time()

  if (getME(smallModel, "is_REML"))
    smallModel <- update(smallModel,REML=FALSE)  
  if (getME(largeModel, "is_REML"))
    largeModel <- update(largeModel,REML=FALSE)
  
  .is.cluster <- !is.null(cl) && inherits(cl, "cluster")
  
  if (!.is.cluster){
    ref <- .getref(largeModel, smallModel, nsim=nsim, seed=seed)
  } else {      
    nsim.cl <- unlist(lapply(splitIndices(nsim, length(cl)), length))
    clusterCall(cl, function() {library(lme4)})
    clusterSetRNGStream(cl)
    ref <- unlist(parallel::clusterMap(cl, .getref, c(largeModel), c(smallModel), nsim.cl))
  }
  
  ctime <- (proc.time()-t0)[3]
  attr(ref,"ctime")   <- ctime
  attr(ref,"samples") <- c(nsim=nsim, npos=sum(ref>0))
  if (details>0)
    cat(sprintf("Reference distribution with %i samples; computing time: %5.2f secs. \n",
                length(ref), ctime))
  
  ref
}



    ##ref <- foo(largeModel, smallModel, cl, nsim)
    ##print(ref)


    ##clusterCall(cl, function() {library(lme4)})
    ##print(clusterCall(cl, function() packageVersion("lme4")))
    ##ref <- unlist(clusterApply(cl, nsim22, function(ii){fun(largeModel, smallModel, nsim=ii)}))
    ##clusterCall(cl, function() detach(package:lme4, unload=TRUE))


##     ## is lme4 loaded:
##     xx <- grep("package:lme4", clusterCall(cl, search))
##     print(xx)
##     ## if yes, then detach
##     if (length(xx)){
##       print(unlist(clusterCall(cl, function() packageVersion("lme4")))); cat("\n")
##       cat("detaching lme4 on clusters\n")
##       clusterCall(cl, function() detach(package:lme4, unload=TRUE))
##     }

##     Sys.sleep(2)

##     cat("load lme4 (the relevant version)\n")
##     clusterCall(cl, function(lll) {library(lme4, lib.loc=lll)}, lib.loc)

## Sys.sleep(2)
##     cat("checking versions:\n")
##     print(unlist(clusterCall(cl, function() packageVersion("lme4")))); cat("\n")
