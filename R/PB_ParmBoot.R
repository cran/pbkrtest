##########################################################
###
### Bartlett corrected LRT
###
##########################################################

PBmodcomp <- function(largeModel, smallModel, nsim=200, ref=NULL, cl=NULL, details=0){
  UseMethod("PBmodcomp")
}

PBmodcomp.mer <- function(largeModel, smallModel, nsim=200, ref=NULL, cl=NULL, details=0){

  if (is.null(ref)){
    ref <- PBrefdist.mer(largeModel, smallModel, nsim=nsim, cl=cl, details=details)
  } else {
    nsim <- length(ref)
  }
  
  lrt     <- .LRT_mer(largeModel, smallModel)
  tobs    <- unname(lrt['tobs'])
  ndf     <- unname(lrt['df'])
  EE      <- mean(ref)
  VV      <- var(ref)

  ##cat(sprintf("EE=%f VV=%f\n", EE, VV))
  
  ## Direct computation of tail probability
  n.extreme <- sum(tobs < ref)
  p.PB  <- n.extreme / length(ref)
  
  ## Bartlett correction - X2 distribution
  BCstat  <- ndf * tobs/EE
  ##cat(sprintf("BCval=%f\n", ndf/EE))
  p.BC    <- 1-pchisq(BCstat,df=ndf)
  
  ## Fit to gamma distribution
  scale   <- VV/EE
  shape   <- EE^2/VV
  p.Ga    <- 1-pgamma(tobs, shape=shape, scale=scale)

  ## Kernel density estimate
  dd <- density(ref)
  p.KD <- sum(dd$y[dd$x>=tobs])/sum(dd$y)

  ## Fit to F-distribution
##   ddf  <- 2*EE/(EE-1)     
##   Fobs <- tobs
##   p.FF <- 1-pf(Fobs, df1=ndf, df2=ddf)

##   rho   <- VV/(2*EE^2)
##   ddf2  <- (ndf*(4*rho+1) - 2)/(rho*ndf-1)
##   lam2  <- (ddf/(ddf-2))/EE
##   Fobs2 <- lam2 * tobs
##   if (ddf2>0)
##     p.FF2 <- 1-pf(Fobs, df1=ndf, df2=ddf2)
##   else
##     p.FF2 <- NA

  # Fit T/d to F-distribution

  ddf  <- 2*EE/(EE-ndf)     
  Fobs <- tobs/ndf
  p.FF <- 1-pf(Fobs, df1=ndf, df2=ddf)
##   rho   <- VV/(2*EE^2)
##   ddf2  <- (ndf*(4*rho+1) - 2)/(rho*ndf-1)
##   lam2  <- (ddf/(ddf-2))/(EE/ndf)
##   print(lam2)
##   Fobs2 <- lam2 * tobs/ndf
##   if (ddf2>0)
##     p.FF2 <- 1-pf(Fobs2, df1=ndf, df2=ddf2)
##   else
##     p.FF2 <- NA


  


  
  f.large <- formula(largeModel)
  attributes(f.large) <- NULL
  f.small <- formula(smallModel)
  attributes(f.small) <- NULL
  
  names(lrt)[1] <- "stat"
  ans <- list(type="X2test",
              f.large=f.large,
              f.small=f.small,
              test = list(
                LRT      = c(c(lrt), ddf=NA),
                PBtest   = c(stat=tobs,    df=NA,   p.value=p.PB,  ddf=NA),
                PBkd     = c(stat=tobs,    df=NA,   p.value=p.KD,  ddf=NA),
                Bartlett = c(stat=BCstat,  df=ndf,  p.value=p.BC,  ddf=NA),
                Gamma    = c(stat=tobs,    df=NA,   p.value=p.Ga,  ddf=NA),
                F        = c(stat=Fobs,    df=ndf,  p.value=p.FF,  ddf=ddf)
#                F2       = c(stat=Fobs2,   df=ndf,  p.value=p.FF2, ddf=ddf2)
                )
              )

  attr(ans,"moment") <- c(mean=EE, var=VV, nsim=nsim)
  attr(ans,"gamma")  <- c(scale=scale, shape=shape)
  attr(ans,"ref")    <- ref
  attr(ans,"ctime")  <- attr(ref,"ctime")
  class(ans) <- c("PBmodcomp", "XXmodcomp")
  ans
}


### ###########################################################
###
### Utilities
###
### ###########################################################

print.XXmodcomp <- function(x, ...){

  cat(sprintf("Parametric bootstrap test; bootstrap samples: %i computing time: %.2f sec.\n",
              attr(x,"moment")['nsim'], round(attr(x,"ctime"),2)  ))
  
  cat("large : "); print(x$f.large)
  cat("small : "); print(x$f.small)
  
                                        #ans      <- as.data.frame(do.call(rbind, x[-c(1:3)]))
  ans      <- as.data.frame(do.call(rbind, x$test))
  ans$p.value    <-round(ans$p.value,options("digits")$digits)
  ans$ddf        <-round(ans$ddf, 3)
  ans$stat       <-round(ans$stat,options("digits")$digits)
  
  print(ans)
}


plot.XXmodcomp <- function(x, ...){

  test <- x$test
  tobs <- test$LRT['stat']
  ref <- attr(x,"ref")
  rr  <- range(ref)
  xx  <- seq(rr[1],rr[2],0.1)
  dd  <- density(ref)
  sc  <- var(ref)/mean(ref)
  sh  <- mean(ref)^2/var(ref)
  
  hist(ref, prob=TRUE,nclass=20, main="Reference distribution")
  abline(v=tobs)
  lines(dd, lty=2, col=2, lwd=2)
  lines(xx,dchisq(xx,df=test$LRT['df']), lty=3, col=3, lwd=2)
  lines(xx,dgamma(xx,scale=sc, shape=sh), lty=4, col=4, lwd=2)
  lines(xx,df(xx,df1=test$F['df'], df2=test$F['ddf']), lty=5, col=5, lwd=2)
  
  smartlegend(x = 'right', y = 'top',
              legend = c("kernel density", "chi-square", "gamma","F"),
              col = 2:5, lty = 2:5)
  
}

as.data.frame.XXmodcomp <- function(x, row.names = NULL, optional = FALSE, ...){
    as.data.frame(do.call(rbind, x[-c(1:3)]))
}

### ###########################################################
###
### Parallel computing of reference distribution
###
### ###########################################################


PBrefdist <- function(largeModel, smallModel, nsim=200, cl=NULL, details=0){
    UseMethod("PBrefdist")
}

PBrefdist.mer <- function(largeModel, smallModel, nsim=200, cl=NULL, details=0){

    t0 <- proc.time()

    if(is.null(smallModel@call$REML) || smallModel@call$REML)
        smallModel <- update(smallModel,REML=FALSE)

    if(is.null(largeModel@call$REML) || largeModel@call$REML)
        largeModel <- update(largeModel,REML=FALSE)

    if (exists(".PBcluster", envir=.GlobalEnv) && inherits(.PBcluster, "cluster")){
      if (details>0)
        cat(sprintf("Using an existing cluster\n"))
      cl <- .PBcluster
    }

    
    .is.null    <- is.null(cl)
    .is.numeric <- is.numeric(cl)
    if (inherits(cl, "cluster")){
      .is.cluster <- TRUE
      .local.cluster <- FALSE
    } else {
      if (is.numeric(cl) && cl > 2){
        cl <- makeSOCKcluster(rep("localhost", round(cl)))
        clusterEvalQ(cl, library(lme4))
        clusterSetupSPRNG(cl)
        .is.cluster <- TRUE
        .local.cluster <- TRUE        
      } else {
        .is.cluster <- FALSE
        .local.cluster <- FALSE
      }
    }
    
    #print(c(.is.cluster=.is.cluster, .local.cluster=.local.cluster))
    
    if (!.is.cluster){
      simdata   <- as.matrix(simulate(smallModel,nsim=nsim))
      ref       <- rep(NA, nsim)
      
      for (kk in seq_len(nsim)){
        yyy      <- simdata[,kk]
        small    <- suppressWarnings(refit(smallModel, newresp=yyy))
        large    <- suppressWarnings(refit(largeModel, newresp=yyy))        
        ttt      <- as.numeric(2*(logLik(large, REML=FALSE) -
                                  logLik(small, REML=FALSE)))
        ref[kk]  <- ttt
      }
    } else {
      
      nsim2 <- round(nsim/length(cl))
      xxx <- clusterCall(cl,
                         function(ll,ss,nsim=200){
                           simdata   <- as.matrix(simulate(ss,nsim=nsim))
                           ref     <- rep(NA, nsim)
                           for (kk in seq_len(nsim)){
                             yyy   <- simdata[,kk]
                             small <- suppressWarnings(refit(ss, newresp=yyy))
                             large <- suppressWarnings(refit(ll, newresp=yyy))
                             ttt   <- as.numeric(2*(logLik(large, REML=FALSE) -
                                                    logLik(small, REML=FALSE)))
                             ref[kk]  <- ttt
                           }
                           ref
                         },
                         largeModel, smallModel, nsim2
                         )
      
      ref <- c(xxx, recursive=TRUE)
      if (.local.cluster)
        stopCluster(cl)
    }
    
    ref <- ref[ref>0]
    ctime <- (proc.time()-t0)[3]
    attr(ref,"ctime") <- ctime
    if (details>0)
      cat(sprintf("Reference distribution with %i samples; computing time: %5.2f secs. \n",
                  length(ref), ctime))
    
    ref
}

makePBcluster <- function(n=2){
  if (exists(".PBcluster", envir=.GlobalEnv) && inherits(.PBcluster, "cluster")){
    cat(sprintf("A cluster exists\n"))
  } else {
    #.PBcluster <<- makeSOCKcluster(rep("localhost", n))
    assign(".PBcluster", makeSOCKcluster(rep("localhost", n)), envir=.GlobalEnv)  
    clusterEvalQ(.PBcluster, library(lme4))
    clusterSetupSPRNG(.PBcluster)
  }
}

stopPBcluster <- function(){
  if (exists(".PBcluster", envir=.GlobalEnv) && inherits(.PBcluster, "cluster")){
    stopCluster(.PBcluster)
    rm(.PBcluster, envir=.GlobalEnv)
  }
}



##########################################################
###
### Likelihood ratio statistic
###
##########################################################

## LRT <- function(largeModel, smallModel){
##   UseMethod("LRT")
## }

.LRT_mer <- function(largeModel, smallModel){
  ll.small <- logLik(smallModel, REML=FALSE)
  ll.large <- logLik(largeModel, REML=FALSE)
  tobs     <- max(0,2*(ll.large-ll.small))
  df11     <- attr(ll.large,"df") - attr(ll.small,"df")
  p.X2     <- 1-pchisq(tobs, df11)
  c(tobs=tobs, df=df11, p.value=p.X2)
}





##########################################################
###
### F-approximated scaled LRT
###
##########################################################

## .FFmodcomp <- function(largeModel, smallModel, nsim=200, ref=NULL, cl=NULL, details=0){
##   UseMethod(".FFmodcomp")
## }

## .FFmodcomp.mer <- function(largeModel, smallModel, nsim=200, ref=NULL, cl=NULL, details=0){

##   if (is.null(ref))
##       ref <- PBrefdist.mer(largeModel, smallModel, nsim=nsim, cl=cl, details=details)

##   lrt   <- .LRT_mer(largeModel, smallModel)
##   ndf <- unname(lrt['df'])
##   names(lrt)[1] <- "tobs"
##   tobs <- unname(lrt['tobs'])

##   m.tref  <- mean(ref)
##   v.tref  <- var(ref)

##   df2        <- round(2 * m.tref/ndf / ( m.tref/ndf - 1 ), 2)
##   df2.ratio  <- df2/(df2-2)
##   FFstat     <- tobs/ndf
##   FFstat2    <- df2/(df2-2)* tobs/(m.tref)

##   p.FF     <- 1-pf(FFstat,  ndf, df2)
##   p.FF2    <- 1-pf(FFstat2, ndf, df2)

##   ans <- list(type="Ftest",
##               f.large=formula(largeModel),
##               f.small=formula(smallModel),
##               LRT = c(lrt,df2=NA),
##               F =c(stat=FFstat,  df=ndf, p=p.FF, df2=df2),
##               F2=c(stat=FFstat2, df=ndf, p=p.FF2, df2=df2)
##               )

##   attr(ans,"moment") <- c(mean=m.tref, var=v.tref)
##   class(ans) <- "XXmodcomp"
##   ans
## }



##########################################################
###
### Parametric bootstrap anova
###
##########################################################

## PBmodcomp <- function(largeModel, smallModel, nsim=200, ref=NULL, cl=NULL, details=0){
##   UseMethod("PBmodcomp")
## }

## PBmodcomp.mer <- function(largeModel, smallModel, nsim=200, ref=NULL, cl=NULL, details=0){

##   if (is.null(ref)){
##     ref <- PBrefdist.mer(largeModel, smallModel, nsim=nsim, cl=cl, details=details)
##   } else {
##     nsim <- length(ref)
##   }
  
##   lrt   <- .LRT_mer(largeModel, smallModel)
##   tobs  <- unname(lrt['tobs'])
##   ndf   <- unname(lrt['df'])

##   m.tref  <- mean(ref)
##   v.tref  <- var(ref)

##   n.extreme <- sum(tobs < ref)
##   p.PB  <- n.extreme / length(ref)

##   f.large <- formula(largeModel)
##   attributes(f.large) <- NULL
##   f.small <- formula(smallModel)
##   attributes(f.small) <- NULL

##   names(lrt)[1] <- "stat"  
##   ans <- list(type="PBtest",  f.large=f.large,   f.small=f.small,
##               LRT = c(lrt),
##               PBtest=c(stat=tobs,    df=NA, p.value=p.PB)
##               )
##   attr(ans,"moment") <- c(mean=m.tref, var=v.tref, nsim=nsim)
##   class(ans) <- "XXmodcomp"
##   ans
## }

