## FRA TIDLIGERE VERSION - og denne virker faktisk!!

## Denne version er der tilføjet diverse i
## vcovAdj2 <- function (object, details = 0)
## {
##     t0 <- proc.time()
##     if (!.is.lmm(object)) {
##         cat("Error in vccovAdj\n")
##         cat(sprintf("model is not a linear mixed moxed model fitted with lmer\n"))
##         stop()
##     }
##     DB <- details > 0
##     if (!(getME(object, "is_REML"))) {
##         object <- update(object, . ~ ., REML = TRUE)
##     }
##     X <- getME(object, "X")
##     Phi <- vcov(object)
##     SigmaG <- LMM_Sigma_G(object, details)
##     ##SG2 <<- SigmaG
##                                         #print(SigmaG)
##     SigmaInv <- chol2inv(chol(forceSymmetric(SigmaG$Sigma)))
##     if (DB) {
##         cat(sprintf("Finding SigmaInv: %10.5f\n", (proc.time() -
##             t0)[1]))
##         t0 <- proc.time()
##     }

##     n.ggamma <- SigmaG$n.ggamma
##     TT <- SigmaInv %*% X
##     HH <- OO <- vector("list", n.ggamma)
##     for (ii in 1:n.ggamma) {
##         .DUM <- SigmaG$G[[ii]] %*% SigmaInv
##         HH[[ii]] <- .DUM
##         OO[[ii]] <- .DUM %*% X
##     }
##     if (DB) {
##         cat(sprintf("Finding TT,HH,OO  %10.5f\n", (proc.time() -
##             t0)[1]))
##         t0 <- proc.time()
##     }
##     ## if(DB){
##     ##     cat("HH:\n"); print(HH); HH2 <<- HH
##     ##     cat("OO:\n"); print(OO); OO2 <<- OO
##     ## }

##     P <- Q <- NULL
##     for (rr in 1:n.ggamma) {
##         OrTrans <- t(OO[[rr]])
##         P <- c(P, list(forceSymmetric(-1 * OrTrans %*% TT)))
##         for (ss in rr:n.ggamma) {
##             Q <- c(Q, list(OrTrans %*% SigmaInv %*% OO[[ss]]))
##         }
##     }
##     if (DB) {
##         cat(sprintf("Finding P,Q:      %10.5f\n", (proc.time() -
##             t0)[1]))
##         t0 <- proc.time()
##     }
##     ## if(DB){
##     ##     cat("P:\n"); print(P); P2 <<- P
##     ##     cat("Q:\n"); print(Q); Q2 <<- Q
##     ## }

##     Ktrace <- matrix(NA, n.ggamma, n.ggamma)
##     for (rr in 1:n.ggamma) {
##         HrTrans <- t(HH[[rr]])
##         for (ss in rr:n.ggamma) {
##             Ktrace[rr, ss] <- Ktrace[ss, rr] <- sum(HrTrans *
##                 HH[[ss]])
##         }
##     }
##     if (DB) {
##         cat(sprintf("Finding Ktrace:   %10.5f\n", (proc.time() -
##             t0)[1]))
##         t0 <- proc.time()
##     }
##     IE2 <- matrix(NA, n.ggamma, n.ggamma)
##     for (ii in 1:n.ggamma) {
##         Phi.P.ii <- Phi %*% P[[ii]]
##         for (jj in c(ii:n.ggamma)) {
##             IE2[ii, jj] <- IE2[jj, ii] <- Ktrace[ii, jj] - 2 *
##                 sum(Phi * Q[[.indexSymmat2vec(ii, jj, n.ggamma)]]) +
##                 sum(Phi.P.ii * (P[[jj]] %*% Phi))
##         }
##     }
##     if (DB) {
##         cat(sprintf("Finding IE2:      %10.5f\n", (proc.time() -
##             t0)[1]))
##         t0 <- proc.time()
##     }
##     eigenIE2 <- eigen(IE2, only.values = TRUE)$values
##     condi <- min(abs(eigenIE2))
##     W <- if (condi > 1e-10)
##         forceSymmetric(2 * solve(IE2))
##     else forceSymmetric(2 * ginv(IE2))

##     ##    print("vcovAdj2")
##     U <- matrix(0, ncol(X), ncol(X))
##     ##print(U)
##     for (ii in 1:(n.ggamma - 1)) {
##         for (jj in c((ii + 1):n.ggamma)) {
##             U <- U + W[ii, jj] * (Q[[.indexSymmat2vec(ii, jj,
##                                                       n.ggamma)]] - P[[ii]] %*% Phi %*% P[[jj]])
##         }
##     }
##     ##print(U)

##     U <- U + t(U)
##     for (ii in 1:n.ggamma) {
##         U <- U + W[ii, ii] * (Q[[.indexSymmat2vec(ii, ii, n.ggamma)]] -
##             P[[ii]] %*% Phi %*% P[[ii]])
##     }
##     ##print(U)

## ##    U2 <<- U
##     GGAMMA <- Phi %*% U %*% Phi
##     PhiA <- Phi + 2 * GGAMMA
##     attr(PhiA, "P") <- P
##     attr(PhiA, "W") <- W
##     attr(PhiA, "condi") <- condi
##     PhiA
## }


## Denne version er helt uberørt

## vcovAdj2 <- function (object, details = 0)
## {
##     if (!.is.lmm(object)) {
##         cat("Error in vccovAdj\n")
##         cat(sprintf("model is not a linear mixed moxed model fitted with lmer\n"))
##         stop()
##     }
##     DB <- details > 0
##     if (!(getME(object, "is_REML"))) {
##         object <- update(object, . ~ ., REML = TRUE)
##     }
##     X <- getME(object, "X")
##     Phi <- vcov(object)
##     SigmaG <- LMM_Sigma_G(object, details)
##     SigmaInv <- chol2inv(chol(forceSymmetric(SigmaG$Sigma)))
##     if (DB) {
##         cat(sprintf("Finding SigmaInv: %10.5f\n", (proc.time() -
##             t0)[1]))
##         t0 <- proc.time()
##     }
##     n.ggamma <- SigmaG$n.ggamma
##     TT <- SigmaInv %*% X
##     HH <- OO <- vector("list", n.ggamma)
##     for (ii in 1:n.ggamma) {
##         .DUM <- SigmaG$G[[ii]] %*% SigmaInv
##         HH[[ii]] <- .DUM
##         OO[[ii]] <- .DUM %*% X
##     }
##     if (DB) {
##         cat(sprintf("Finding TT,HH,OO  %10.5f\n", (proc.time() -
##             t0)[1]))
##         t0 <- proc.time()
##     }
##     P <- Q <- NULL
##     for (rr in 1:n.ggamma) {
##         OrTrans <- t(OO[[rr]])
##         P <- c(P, list(forceSymmetric(-1 * OrTrans %*% TT)))
##         for (ss in rr:n.ggamma) {
##             Q <- c(Q, list(OrTrans %*% SigmaInv %*% OO[[ss]]))
##         }
##     }
##     if (DB) {
##         cat(sprintf("Finding P,Q:      %10.5f\n", (proc.time() -
##             t0)[1]))
##         t0 <- proc.time()
##     }
##     Ktrace <- matrix(NA, n.ggamma, n.ggamma)
##     for (rr in 1:n.ggamma) {
##         HrTrans <- t(HH[[rr]])
##         for (ss in rr:n.ggamma) {
##             Ktrace[rr, ss] <- Ktrace[ss, rr] <- sum(HrTrans *
##                 HH[[ss]])
##         }
##     }
##     if (DB) {
##         cat(sprintf("Finding Ktrace:   %10.5f\n", (proc.time() -
##             t0)[1]))
##         t0 <- proc.time()
##     }
##     IE2 <- matrix(NA, n.ggamma, n.ggamma)
##     for (ii in 1:n.ggamma) {
##         Phi.P.ii <- Phi %*% P[[ii]]
##         for (jj in c(ii:n.ggamma)) {
##             IE2[ii, jj] <- IE2[jj, ii] <- Ktrace[ii, jj] - 2 *
##                 sum(Phi * Q[[.indexSymmat2vec(ii, jj, n.ggamma)]]) +
##                 sum(Phi.P.ii * (P[[jj]] %*% Phi))
##         }
##     }
##     if (DB) {
##         cat(sprintf("Finding IE2:      %10.5f\n", (proc.time() -
##             t0)[1]))
##         t0 <- proc.time()
##     }
##     eigenIE2 <- eigen(IE2, only.values = TRUE)$values
##     condi <- min(abs(eigenIE2))
##     W <- if (condi > 1e-10)
##         forceSymmetric(2 * solve(IE2))
##     else forceSymmetric(2 * ginv(IE2))
##     U <- matrix(0, ncol(X), ncol(X))
##     for (ii in 1:(n.ggamma - 1)) {
##         for (jj in c((ii + 1):n.ggamma)) {
##             U <- U + W[ii, jj] * (Q[[.indexSymmat2vec(ii, jj,
##                 n.ggamma)]] - P[[ii]] %*% Phi %*% P[[jj]])
##         }
##     }
##     U <- U + t(U)
##     for (ii in 1:n.ggamma) {
##         U <- U + W[ii, ii] * (Q[[.indexSymmat2vec(ii, ii, n.ggamma)]] -
##             P[[ii]] %*% Phi %*% P[[ii]])
##     }
##     GGAMMA <- Phi %*% U %*% Phi
##     PhiA <- Phi + 2 * GGAMMA
##     attr(PhiA, "P") <- P
##     attr(PhiA, "W") <- W
##     attr(PhiA, "condi") <- condi
##     PhiA
## }
