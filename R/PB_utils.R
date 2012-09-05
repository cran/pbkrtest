
.formula2list <- function(form){
  lhs <- form[[2]]  
  tt  <- terms(form)
  tl  <- attr(tt, "term.labels")
  r.idx <- grep("\\|", tl)

  if (length(r.idx)){
    rane  <- paste("(", tl[r.idx], ")")
    f.idx <- (1:length(tl))[-r.idx]
    if (length(f.idx))
      fixe  <- tl[f.idx]
    else
      fixe  <- NULL
  } else {
    rane <- NULL 
    fixe <- tl
  }
  
  ans <- list(lhs=deparse(lhs),
              rhs.fix=fixe,
              rhs.ran=rane)
  ans
}


restrictionMatrix2model <- function(largeModel, LL){
  UseMethod("restrictionMatrix2model")
}

restrictionMatrix2model.mer <- function(largeModel, LL){

  form <- as.formula(formula(largeModel))
  XX.lg 	 <- getME(largeModel, "X")
  
  attributes(XX.lg)[-1] <- NULL
  XX.sm <- .restrictedModelMatrix(XX.lg, LL)
  
  ncX.sm  <- ncol(XX.sm)
  colnames(XX.sm) <- paste(".X", 1:ncX.sm, sep='')
  
  rhs.fix2 <- paste(".X", 1:ncX.sm, sep='', collapse="+")
  
  fff  <- .formula2list(form)
  new.formula <- as.formula(paste(fff$lhs, "~ -1+", rhs.fix2, "+", fff$rhs.ran))
  new.data    <- cbind(XX.sm, eval(largeModel@call$data))
  
##   ans <- lmer(eval(new.formula), data=new.data, REML=getME(largeModel, "is_REML"))
  ans <- update(largeModel, eval(new.formula), data=new.data)
  ans
}

restrictionMatrix2model.lm <- function(largeModel, LL){

  form <- as.formula(formula(largeModel))
  XX.lg 	 <- model.matrix(largeModel)
  attributes(XX.lg)[-1] <- NULL
  XX.sm <- .restrictedModelMatrix(XX.lg, LL)
  
  ncX.sm  <- ncol(XX.sm)
  colnames(XX.sm) <- paste(".X", 1:ncX.sm, sep='')
  
  rhs.fix2 <- paste(".X", 1:ncX.sm, sep='', collapse="+")  
  fff  <- .formula2list(form)
  new.formula <- as.formula(paste(fff$lhs, "~ -1+", rhs.fix2))
  new.data    <- cbind(XX.sm, eval(largeModel$data))
  #print(new.data)
  ans <- update(largeModel, eval(new.formula), data=new.data)
  ans
}



##########################################################
###
### Likelihood ratio statistic
###
##########################################################

getLRT <- function(largeModel, smallModel){
  UseMethod("getLRT")
}

getLRT.mer <- function(largeModel, smallModel){
  ll.small <- logLik(smallModel, REML=FALSE)
  ll.large <- logLik(largeModel, REML=FALSE)
  tobs     <- 2*(ll.large-ll.small)
  df11     <- attr(ll.large, "df") - attr(ll.small, "df")
  p.X2     <- 1-pchisq(tobs, df11)
  c(tobs=tobs, df=df11, p.value=p.X2)
}


getLRT.lm <- function(largeModel, smallModel){
  ll.small <- logLik(smallModel)
  ll.large <- logLik(largeModel)
  tobs     <- 2*(ll.large-ll.small)
  df11     <- attr(ll.large, "df") - attr(ll.small, "df")
  p.X2     <- 1-pchisq(tobs, df11)
  c(tobs=tobs, df=df11, p.value=p.X2)
}







