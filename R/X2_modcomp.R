## ##########################################################################
##
#' @title Chisq test 
#' 
#' @description Chisq test 
#' @concept model_comparison
#' @name x2__modcomp
#' 
## ##########################################################################
#' @details TBW
#' 
#' @param largeModel An \code{lmer} model
#' @param smallModel An \code{lmer} model or a restriction matrix
#' @param betaH A number or a vector of the beta of the hypothesis,
#'     e.g. L beta=L betaH. If `smallModel` is a model object then betaH=0.
#' @param details If larger than 0 some timing details are printed.
#' @param ... Additional arguments, currently not used.
#' 
#' @author Ulrich Halekoh \email{uhalekoh@@health.sdu.dk}, Søren Højsgaard
#'     \email{sorenh@@math.aau.dk}
#'
#' (fm0 <- lmer(Reaction ~ (Days|Subject), sleepstudy))
#' (fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
#' (fm2 <- lmer(Reaction ~ Days + I(Days^2) + (Days|Subject), sleepstudy))
#'
#' ## Test for no effect of Days in fm1, i.e. test fm0 under fm1
#' 
#' X2modcomp(fm1, "Days")
#' X2modcomp(fm1, ~.-Days)
#' L1 <- cbind(0, 1)
#' ## X2modcomp(fm1, L1) ## FIXME

#' X2modcomp(fm1, fm0)
#' anova(fm1, fm0)


#' @export
#' @rdname x2__modcomp
X2modcomp <- function(largeModel, smallModel, betaH=0, details=0, ...){
    UseMethod("X2modcomp")
}

#' @export
#' @rdname x2__modcomp
X2modcomp.default <- function(largeModel, smallModel, betaH=0, details=0, ...) {
    X2modcomp_internal(largeModel=largeModel, smallModel=smallModel, betaH=betaH, details=details)
}


X2modcomp_internal <- function(largeModel, smallModel, betaH=0, details=0) {

    if (is.character(smallModel))
        smallModel <- doBy::formula_add_str(formula(largeModel), terms=smallModel, op="-")
    
    if (inherits(smallModel, "formula"))
        smallModel  <- update(largeModel, smallModel)

    ## if (is.numeric(smallModel) && !is.matrix(smallModel))
        ## smallModel <- matrix(smallModel, nrow=1)
    
    ## w <- modcomp_init(largeModel, smallModel, matrixOK = TRUE)

    ## if (w == -1) stop('Models have equal mean stucture or are not nested!')
    ## if (w == 0){
    ##     ## First given model is submodel of second; exchange the models
    ##     tmp <- largeModel; 
    ##     largeModel <- smallModel; 
    ##     smallModel <- tmp
    ## }
    
    ## ## Refit large model with REML if necessary
    ## if (!(getME(largeModel, "is_REML"))){
    ##     largeModel <- update(largeModel, .~., REML=TRUE)
    ## }

    ## print(largeModel)
    ## print(smallModel)
    
    X2modcomp_worker(largeModel, smallModel, betaH=betaH, details=details)    
}



X2modcomp_worker <- function(largeModel, smallModel, betaH=0, details=0) {

    if (is.null(betaH)) betaH <- 0
    if (is.null(details)) details <- 0
    
    ## cat("X2modcomp_worker\n")
    ## All computations are based on 'largeModel' and the restriction matrix 'L'
    ## -------------------------------------------------------------------------
    t0    <- proc.time()
    L     <- NULL ##model2restriction_matrix(largeModel, smallModel)

    ## if (inherits(smallModel, "matrix")){
        ## smallModel <- suppressWarnings(restriction_matrix2model(largeModel, L=smallModel))
    ## }


    ## PhiA  <- vcovAdj(largeModel, details)
    ## stats <- .KR_adjust(PhiA, Phi=vcov(largeModel), L, beta=fixef(largeModel), betaH)
    ## stats <- lapply(stats, c) ## To get rid of all sorts of attributes

    stats <- NULL

    LRTstat     <- getLRT(largeModel, smallModel)  
    ## cat("LRTstat:\n"); print(LRTstat); cat("LRTstat done:\n");
    
    ans   <- X2compute_p_values(LRTstat,stats)
    ## print(ans)
    
    formula.large <- formula(largeModel)
    formula.small <- formula(smallModel)    
    attributes(formula.large) <- NULL
    
    ans$formula.large <- formula.large
    ans$formula.small <- formula.small
    ans$ctime   <- (proc.time() - t0)[3]
    ans$L       <- L

    ## print(ans)
    
    out <- ans$test[1,, drop=FALSE]
    ## print(out)
    attr(out, "aux") <- ans

    attr(out, "heading") <- c(
        deparse(formula.large),
        deparse(formula.small))

    class(out) <- c("X2modcomp", "anova", "data.frame")
    return(out)
}

X2compute_p_values <- function(LRTstat, stats=NULL){

    tobs <- unname(LRTstat[1])
    ndf  <- unname(LRTstat[2])
    p.chi <- 1 - pchisq(tobs, df=ndf)
    
    test = list(
        LRT        = c(stat=tobs,            df=ndf,        ddf=NA,                       p.value=p.chi)
        )
    test  <- as.data.frame(do.call(rbind, test))
    test$df <- as.integer(test$df)
    out   <- list(test=test, type="X2", aux=stats$aux, stats=stats)
    ## Notice: stats are carried to the output. They are used for get getKR function...
    class(out) <- c("X2modcomp")
    out
}


## #' @export
## #' @rdname x2__modcomp
## X2modcomp.lmerMod <- function(largeModel, smallModel, betaH=0, details=0, ...) {
##     X2modcomp_internal(largeModel=largeModel, smallModel=smallModel, betaH=betaH, details=details)
## }

## #' @export
## #' @rdname x2__modcomp
## X2modcomp.glmerMod <- function(largeModel, smallModel, betaH=0, details=0, ...) {
##     X2modcomp_internal(largeModel=largeModel, smallModel=smallModel, betaH=betaH, details=details)
## }

## #' @export
## #' @rdname x2__modcomp
## X2modcomp.gls <- function(largeModel, smallModel, betaH=0, details=0, ...) {
##     X2modcomp_internal(largeModel=largeModel, smallModel=smallModel, betaH=betaH, details=details)
## }

## #' @export
## #' @rdname x2__modcomp
## X2modcomp.lm <- function(largeModel, smallModel, betaH=0, details=0, ...) {
##     X2modcomp_internal(largeModel=largeModel, smallModel=smallModel, betaH=betaH, details=details)
## }
