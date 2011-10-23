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
  m.tref  <- mean(ref)
  v.tref  <- var(ref)


  ## Direct computation of tail probability
  n.extreme <- sum(tobs < ref)
  p.PB  <- n.extreme / length(ref)
  
  ## Bartlett correction - X2 distribution
  BCstat  <- ndf * tobs/m.tref
  p.BC    <- 1-pchisq(BCstat,df=ndf)
  
  ## Fit to gamma distribution
  scale   <- v.tref/m.tref
  shape   <- m.tref^2/v.tref
  p.Ga    <- 1-pgamma(tobs, shape=shape, scale=scale)
  
  res <- list(c(lrt, p.BC=p.BC, BCstat=BCstat, m.tref=m.tref),
              f.large=formula(largeModel), f.small=formula(smallModel))

  f.large <- formula(largeModel)
  attributes(f.large) <- NULL
  f.small <- formula(smallModel)
  attributes(f.small) <- NULL
  
  names(lrt)[1] <- "stat"
  ans <- list(type="X2test", f.large=f.large, f.small=f.small,
              LRT      = c(lrt),
              PBtest   = c(stat=tobs,    df=NA,  p.value=p.PB),
              Bartlett = c(stat=BCstat,  df=ndf, p.value=p.BC),
              Gamma    = c(stat=tobs,    df=NA,  p.value=p.Ga)
              )
  attr(ans,"moment") <- c(mean=m.tref, var=v.tref, nsim=nsim)
  attr(ans,"gamma")  <- c(scale=scale, shape=shape)
  attr(ans,"ref")    <- ref
  ##attr(ans,"gamma2") <- c(scale=scale2, shape=shape2)
  class(ans) <- c("PBmodcomp", "XXmodcomp")
  ans
}


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



### ###########################################################
###
### Utilities
###
### ###########################################################


print.XXmodcomp <- function(x, ...){

    cat("large : "); print(x$f.large)
    cat("small : "); print(x$f.small)

    cat(sprintf("Number of parametric bootstrap samples: %i\n",
                attr(x,"moment")['nsim']))
    ans      <- as.data.frame(do.call(rbind, x[-c(1:3)]))
    ans$p.value    <-round(ans$p.value,options("digits")$digits)
    ans$stat <-round(ans$stat,options("digits")$digits)

    print(ans)
}


plot.XXmodcomp <- function(x, ...){

  ref <- attr(x,"ref")
  rr  <- range(ref)
  xx  <- seq(rr[1],rr[2],0.1)
  dd  <- density(ref)
  sc  <- var(ref)/mean(ref)
  sh  <- mean(ref)^2/var(ref)
  
  hist(ref, prob=TRUE,nclass=20, main="Reference distribution")
  abline(v=x$LRT['stat'])
  lines(dd, lty=2, col=2, lwd=2)
  lines(xx,dchisq(xx,df=x$LRT['df']), lty=3, col=3, lwd=2)
  lines(xx,dgamma(xx,scale=sc, shape=sh ), lty=4, col=4, lwd=2)
  
  smartlegend(x = 'center', y = 'top',
              legend = c("kernel density", "chi-square", "gamma"), col = 2:4, lty = 2:4)
  
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
    
    if (details>0)
      cat(sprintf("Reference distribution with %i samples; computing time: %5.2f secs. \n",
                  length(ref), (proc.time()-t0)[3]))
    
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
