##########################################################
###
### Bartlett corrected LRT
###
##########################################################

PBmodcomp <- function(largeModel, smallModel, nsim=200, ref=NULL, cl=NULL, details=0){
  UseMethod("PBmodcomp")
}

PBmodcomp.mer <- function(largeModel, smallModel, nsim=200, ref=NULL, cl=NULL, details=0){

  f.large <- formula(largeModel)
  attributes(f.large) <- NULL

  if (inherits(smallModel, c("Matrix", "matrix"))){
    f.small <- smallModel
    smallModel <- restrictionMatrix2model(largeModel, smallModel)
  } else {
    f.small <- formula(smallModel)
    attributes(f.small) <- NULL
  }
  
  if (is.null(ref)){
    ref <- PBrefdist(largeModel, smallModel, nsim=nsim, cl=cl, details=details)
  } else {
    nsim <- length(ref)
  }
  
  LRTstat     <- getLRT(largeModel, smallModel)
  ans         <- .finalizePB(LRTstat, ref)
  ans$LRTstat <- LRTstat
  ans$ref     <- ref
  ans$f.large <- f.large
  ans$f.small <- f.small
  ans
}

PBmodcomp.lm <- function(largeModel, smallModel, nsim=200, ref=NULL, cl=NULL, details=0){

  ok.fam <- c("binomial","gaussian","Gamma","inverse.gaussian","poisson")
  f.large <- formula(largeModel)
  attributes(f.large) <- NULL

  if (inherits(smallModel, c("Matrix", "matrix"))){
    #cat(".....smallModel is a matrix\n")
    f.small <- smallModel
    smallModel <- restrictionMatrix2model(largeModel, smallModel)
  } else {
    #cat(".....smallModel is a model object\n")
    f.small <- formula(smallModel)
    attributes(f.small) <- NULL
  }

  if (!all.equal((fam.l <- family(largeModel)), (fam.s <- family(smallModel))))
    stop("Models do not have identical identical family\n")
  if (!(fam.l$family %in% ok.fam)){
    stop(sprintf("family must be of type %s", toString(ok.fam)))
  }
  
  if (is.null(ref)){
    ref <- PBrefdist(largeModel, smallModel, nsim=nsim, cl=cl, details=details)
  } else {
    nsim <- length(ref)
  }
  
  LRTstat     <- getLRT(largeModel, smallModel)
  ans         <- .finalizePB(LRTstat,ref)
  ans$LRTstat <- LRTstat
  ans$ref     <- ref
  ans$f.large <- f.large
  ans$f.small <- f.small
  ans
}



.finalizePB <- function(LRTstat, ref){

  tobs <- unname(LRTstat[1])
  ndf  <- unname(LRTstat[2])
  
  nsim    <- length(ref)
  ##cat(sprintf("EE=%f VV=%f\n", EE, VV))
  p.chi <- 1-pchisq(tobs, df=ndf)
  ## Direct computation of tail probability
  n.extreme <- sum(tobs < ref)
  p.PB  <- n.extreme / nsim
  
  test = list(
    LRT      = c(stat=tobs,    df=ndf,    p.value=p.chi),
    PBtest   = c(stat=tobs,    df=NA,     p.value=p.PB))

  test  <- as.data.frame(do.call(rbind, test))
  ans   <- list(test=test, type="X2test", ctime=attr(ref,"ctime"))
  class(ans) <- c("PBmodcomp")
  ans
}

.summarizePB <- function(LRTstat, ref){

  tobs <- unname(LRTstat[1])
  ndf  <- unname(LRTstat[2])
  
  EE      <- mean(ref)
  VV      <- var(ref)
  nsim    <- length(ref)
  ##cat(sprintf("EE=%f VV=%f\n", EE, VV))

  p.chi <- 1-pchisq(tobs, df=ndf)
  
  ## Direct computation of tail probability
  n.extreme <- sum(tobs < ref)
  p.PB  <- n.extreme / nsim

  se <- round(sqrt(p.PB*(1-p.PB)/nsim),4)
  ci <- round(c(-1.96, 1.96)*se + p.PB,4)

  ## Kernel density estimate
  dd <- density(ref)
  p.KD <- sum(dd$y[dd$x>=tobs])/sum(dd$y)
  
  ## Bartlett correction - X2 distribution
  BCstat  <- ndf * tobs/EE
  ##cat(sprintf("BCval=%f\n", ndf/EE))
  p.BC    <- 1-pchisq(BCstat,df=ndf)
  
  ## Fit to gamma distribution
  scale   <- VV/EE
  shape   <- EE^2/VV
  p.Ga    <- 1-pgamma(tobs, shape=shape, scale=scale)

  ## Fit T/d to F-distribution
  ddf  <- 2*EE/(EE-ndf)     
  Fobs <- tobs/ndf
  p.FF <- 1-pf(Fobs, df1=ndf, df2=ddf)

  
  test = list(
    LRT      = c(stat=tobs,    df=ndf, ddf=NA,   p.value=p.chi),
    PBtest   = c(stat=tobs,    df=NA,  ddf=NA,   p.value=p.PB),
    PBkd     = c(stat=tobs,    df=NA,  ddf=NA,   p.value=p.KD),
    Gamma    = c(stat=tobs,    df=NA,  ddf=NA,   p.value=p.Ga),
    F        = c(stat=Fobs,    df=ndf, ddf=ddf,  p.value=p.FF),
    Bartlett = c(stat=BCstat,  df=ndf, ddf=NA,   p.value=p.BC)
    ##F2       = c(stat=Fobs2,   df=ndf,  p.value=p.FF2, ddf=ddf2)
    )

  test <- as.data.frame(do.call(rbind, test))
  ans <- list(test=test, type="X2test")
  
  attr(ans,"moment") <- c(mean=EE, var=VV, nsim=nsim)
  attr(ans,"gamma")  <- c(scale=scale, shape=shape)
  attr(ans,"ref")    <- ref
  attr(ans,"ci")     <- ci
  attr(ans,"se")     <- se
  attr(ans,"ctime")  <- attr(ref,"ctime")
  class(ans) <- c("PBmodcomp")
  ans
}


### ###########################################################
###
### Utilities
###
### ###########################################################

.PBcommon <- function(x){

  cat(sprintf("Parametric bootstrap test; "))
  
  if (!is.null((zz<- attr(x,"moment")))){
    cat(sprintf("bootstrap samples: %i", zz['nsim']))
  }
  if (!is.null((zz<- attr(x,"ctime")))){
    cat(sprintf(" computing time: %.2f sec.", round(zz,2)))
  }
  cat("\n")
  
  if(!is.null(x$f.large)){
    cat("large : "); print(x$f.large)
    cat("small : "); print(x$f.small)
  }
}


print.PBmodcomp <- function(x, ...){
  .PBcommon(x)
  tab <- x$test
  printCoefmat(tab, tst.ind=1, na.print='', has.Pvalue=TRUE)
  return(invisible(x))
}


summary.PBmodcomp <- function(object,...){
  ans <- .summarizePB(object$LRTstat, object$ref)
  ans$f.large <- object$f.large
  ans$f.small <- object$f.small
  class(ans) <- "summaryPB"
  ans
}

print.summaryPB <- function(x,...){
  .PBcommon(x)
  ans <- x$test
  printCoefmat(ans, tst.ind=1, na.print='', has.Pvalue=TRUE)
  cat("\n")

  ci <- attr(x, "ci")  
  #ci <- x$ci
  cat(sprintf("95 pct CI for PBtest   : [%s]\n", toString(ci)))

  mo <- attr(x, "moment") 
  #mo <- x$moment
  cat(sprintf("Reference distribution : mean=%f var=%f nsim=%d\n", mo[1], mo[2], mo[3]))

  ga <- attr(x, "gamma")
  #ga <- x$gamma
  cat(sprintf("Gamma approximation    : scale=%f shape=%f\n", ga[1], ga[2]))

  return(invisible(x))
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






##   rho   <- VV/(2*EE^2)
##   ddf2  <- (ndf*(4*rho+1) - 2)/(rho*ndf-1)
##   lam2  <- (ddf/(ddf-2))/(EE/ndf)
##   cat(sprintf("EE=%f, VV=%f, rho=%f, lam2=%f\n",
##               EE, VV, rho, lam2))

##   ddf2 <- 4 + (ndf+2)/(rho*ndf-1)
  
##   Fobs2 <- lam2 * tobs/ndf
##   if (ddf2>0)
##     p.FF2 <- 1-pf(Fobs2, df1=ndf, df2=ddf2)
##   else
##     p.FF2 <- NA
