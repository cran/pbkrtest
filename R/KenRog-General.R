
KRmodcomp <- function(largeModel, smallModel,beta0=0, details=0){
    UseMethod("KRmodcomp")
}

KRmodcomp.mer<-function(largeModel,smallModel,beta0=0, details=0) {
    ##smallModel can either be a lmer model or a restriction matrix L

    w <- modcomp_init(largeModel,smallModel,matrixOK=TRUE)
    if (w==-1) {
        print ('Error in KRmodcomp')
        print( 'both models have either equal fixed mean stucture')
        print( 'or are not nested')
        stop()
    }
    if (w==0){
        print ('Error in KRmodcomp')
        print( 'first given model is submodel of second')
        print(' exchange the models')
        stop()
    }



    ## refitting large model with REML if necessary
    ## in modcomp_inmitmit is checked that legreModel is  a gaussian gaussian mixed model 
    largeModel<-
        if (largeModel@dims['REML'] == 1)
        {
            largeModel
        }
        else
        {
		    warning("\n largeModel has been refitted with REML=TRUE \n")
            update(largeModel,.~.,REML=TRUE)
        }


    L<- .createRestrictionMatrix(largeModel,smallModel)

    
    ## All further computations are based on 'largeModel' and the restriction matrix 'L'
    ## -------------------------------------------------------------------------

    t0 <- proc.time()
    stats <- .KRmodcompPrimitive(largeModel, L, beta0, details)
    formSmall <-
      if ('mer' %in% class(smallModel)){
        .zzz <- formula(smallModel)
        attributes(.zzz) <- NULL
        .zzz
      } else {
        list(L=L,beta0=beta0)
      }
    formLarge <- formula(largeModel)
    attributes(formLarge) <- NULL
    
    res<-list(stats=stats,f.large=formLarge,f.small=formSmall)
    attr(res,"ctime")  <- (proc.time()-t0)[1]
    class(res)<-c("KRmodcomp")
    res
}




.KRmodcompPrimitive<-function(largeModel, L, beta0, details) {
  PhiA<-vcovAdj(largeModel, details)
  .KR_adjust(PhiA, Phi=vcov(largeModel), L, beta=fixef(largeModel), beta0 )
}

vcovAdj <- function(largeModel, details=0) {
  if(!.is.lmm(largeModel)) {
    cat("Error in vccovAdj\n")
    cat(sprintf("largeModel is not a linear mixed moxed model fittedt with lmer\n"))
    stop()
  }
  DB <- details>0

    ## refitting large model with REML if necessary
    ## in modcomp_inmitmit is checked that legreModel is  a gaussian gaussian mixed model 
    largeModel<-
        if (largeModel@dims['REML'] == 1)
        {
            largeModel
        }
        else
        {
		    warning("\n largeModel has been refitted with REML=TRUE \n")
            update(largeModel,.~.,REML=TRUE)
		}

  
  X<-getME(largeModel,"X")
  
  Phi    <- vcov(largeModel)
  GGamma <- VarCorr(largeModel)
                                        # s -> n.varcomp
  n.groupFac<- largeModel@dims['nt'] #= number of random effects terms (..|..)
                                        # (..|F1) + (..|F1) are group factors!
                                        # without the residual variance
  
  ## size of the symmetric variance Gamma_i for reach groupFac
  nn.GGamma <- integer(n.groupFac)
  ggamma <- NULL
  for (ii in 1: (n.groupFac)) {
    Lii<-GGamma[[ii]]
    nu<-ncol(Lii)
    nn.GGamma[ii]<- nu
    ## The lower.tri construxtion esnures, that (because Lii is symmetric!)
    ## Lii[lower.tri(Lii,diag=TRUE)= Lii[1,1],Lii[1,2],Lii[1,3]..Lii[1,nu],
    ##                               Lii[2,2], Lii[2,3] ...
    ggamma<-c(ggamma,Lii[lower.tri(Lii,diag=TRUE)])
  }
  
  ## number of variance parameters of each GGamma_i
  mm.GGamma<- nn.GGamma * (nn.GGamma+1)/2

  ##adding the residuals variance to ggamma
  ##so in ggamma nd n.ggamma the residual variance is included!
  ggamma<-c(ggamma,attr(GGamma,'sc')^2)
  n.ggamma<-length(ggamma)
  
  ##
  group.index<-getME(largeModel,"Gp")
  nn.groupFac<-diff(group.index)
  
  ## number of random effects in each groupFac
  ## residual error here excluded!
  nn.groupFacLevels<-nn.groupFac/nn.GGamma
  
  Zt<-getME(largeModel,"Zt")
  ## G_r:
  
  t0 <- proc.time()
  
  G<-NULL
  for (ss in 1:n.groupFac)
    {
      zIndex.sub<-group.index[ss]+
        1+c(0:(nn.GGamma[ss]-1))*nn.groupFacLevels[ss] +
          rep(0:(nn.groupFacLevels[ss]-1),each=nn.GGamma[ss])
      ##ZZ<-Zt[ (index.nn.group[ss]+1):index.nn.group[ss+1], ]
      ZZ<-Zt[zIndex.sub, ]
                                        #cat("dim(ZZ)"); print(dim(ZZ))
      Ig<-sparseMatrix(1:nn.groupFacLevels[ss],
                       1:nn.groupFacLevels[ss],x=1)
      
      for (rr in 1:mm.GGamma[ss] )
        {
          ii.jj <- .indexVec2Symmat(rr,nn.GGamma[ss])
          ii.jj <- unique(ii.jj)
          EE    <-
            if (length(ii.jj)==1){
              sparseMatrix(ii.jj,ii.jj,x=1,dims=rep(nn.GGamma[ss],2))
            } else {
              sparseMatrix(ii.jj,ii.jj[2:1],dims=rep(nn.GGamma[ss],2))
            }          
          EE <- Ig %x% EE  ## Kronecker product
          G  <- c(G,list(t(ZZ)%*% EE %*% ZZ))
        }
    }
  
  G<-c(G,list(sparseMatrix(1:nrow(X),1:nrow(X),x=1))) ## The last one is for the residual!
  if(DB){cat(sprintf("Finding G         %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}
  
  Sigma<-ggamma[1]*G[[1]]
  for (ii in 2:n.ggamma) {
    Sigma<- Sigma + ggamma[ii] * G[[ii]]
  }
  if(DB){cat(sprintf("Finding Sigma:    %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}
  
  SigmaInv <- chol2inv(chol(forceSymmetric(Sigma)))
  if(DB){cat(sprintf("Finding SigmaInv: %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}
  
  ## Finding, TT, HH, 00
  TT <- SigmaInv %*% X
  HH <- OO <-  vector("list", n.ggamma)
  for (ii in 1:n.ggamma) {
    .DUM<-G[[ii]] %*% SigmaInv
    HH[[ii]] <- .DUM
    OO[[ii]] <- .DUM %*% X
  }
  if(DB){cat(sprintf("Finding TT,HH,OO  %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}
  
  P <- Q <-NULL
  for (rr in 1:n.ggamma) {
    OrTrans <- t(OO[[rr]])
    P <- c(P, list(forceSymmetric( -1 * OrTrans %*%  TT)))
    for (ss in rr:n.ggamma) {
      Q <- c(Q,list(OrTrans %*% SigmaInv %*% OO[[ss]] ))
    }}
  if(DB){cat(sprintf("Finding P,Q:      %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}
  
  Ktrace <- matrix(NA,n.ggamma,n.ggamma)
  for (rr in 1:n.ggamma) {
    HrTrans<-t(HH[[rr]])
    for (ss in rr:n.ggamma) {
      Ktrace[rr,ss] <- Ktrace[ss,rr]<- sum( HrTrans * HH[[ss]])
    }}
  if(DB){cat(sprintf("Finding Ktrace:   %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}
  
  ## Nov. 24. 2011; SHD
  ## Alternative computation of Ktrace. Seems to be no faster than the one above but please do
  ## not delete
  ##     Ktrace2 <- matrix(NA,n.ggamma,n.ggamma)
  ##     for (rr in 1:n.ggamma) {
  ##       HrTrans<-t(H[[rr]])
  ##       Ktrace2[rr,rr] <- sum( HrTrans * t(HrTrans))
  ##       if (rr < n.ggamma){
  ##         for (ss in (rr+1):n.ggamma) {
  ##           Ktrace2[rr,ss] <- Ktrace2[ss,rr]<- sum( HrTrans * H[[ss]])
  ##         }}}
  ##     cat(sprintf("Finding Ktrace(2): %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()
  ##     print(Ktrace)
  ##     print(Ktrace2)
  ##     print(Ktrace2-Ktrace)
  
  IE2 <- matrix(NA,n.ggamma,n.ggamma)
  for (ii in 1:n.ggamma) {
    Phi.P.ii <- Phi %*% P[[ii]]
    for (jj in c(ii:n.ggamma)) {
      IE2[ii,jj]<- IE2[jj,ii]<- Ktrace[ii,jj] - 
        2 * sum(Phi*Q[[.indexSymmat2vec(ii,jj,n.ggamma)]]) +
          sum( Phi.P.ii * (  P[[jj]] %*% Phi))
    }}
  if(DB){cat(sprintf("Finding IE2:      %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}
  
  eigenIE2 <- eigen(IE2,only.values=TRUE)$values
  condi    <- min(abs(eigenIE2))
  
  W <- if(condi>1e-10) forceSymmetric(2* solve(IE2)) else forceSymmetric(2* ginv(IE2))
  
  U<-matrix(0,ncol(X),ncol(X))
  for (ii in 1:(n.ggamma-1)) {
    for (jj in c((ii+1):n.ggamma)) {
      U <- U + W[ii,jj] * (Q[[.indexSymmat2vec(ii,jj,n.ggamma)]] -
                           P[[ii]] %*% Phi %*% P[[jj]])
    }}
  
### FIXME: Ulrich: Er det ikke sådan, at du her får beregnet diagonalen med to gange???
  U<-U+t(U) 
  
  for (ii in 1:n.ggamma) {
    U<- U +   W[ii,ii] * (Q[[.indexSymmat2vec(ii,ii,n.ggamma)]]- P[[ii]] %*% Phi %*% P[[ii]])
  }
  
  GGAMMA <-  Phi %*% U %*% Phi
  PhiA   <-  Phi + 2* GGAMMA
  attr(PhiA,"P")<-P
  attr(PhiA,"W")<-W
  PhiA
}

.KR_adjust <- function(PhiA,Phi,L,beta,beta0){  
  Theta  <-  t(L) %*% solve( L %*% Phi %*% t(L), L)
  P<-attr(PhiA,"P")
  W<-attr(PhiA,"W")
  
  A1<-A2<-0
  ThetaPhi<-Theta%*%Phi
  n.ggamma<-length(P)
  for (ii in 1:n.ggamma) {
    for (jj in c(ii:n.ggamma)) {
      e<-ifelse(ii==jj, 1, 2)
      ui<-ThetaPhi %*% P[[ii]] %*% Phi
      uj<-ThetaPhi %*% P[[jj]] %*% Phi
      A1<- A1+  e* W[ii,jj] * (.spur(ui) * .spur(uj))
      A2<- A2+  e* W[ii,jj] *  sum(ui * t(uj))
    }}
  

  q <- rankMatrix(L)
  B <- (1/(2*q)) * (A1+6*A2)
  g <- ( (q+1)*A1 - (q+4)*A2 )  / ((q+2)*A2)
  c1<- g/(3*q+ 2*(1-g))
  c2<- (q-g) / (3*q + 2*(1-g))
  c3<- (q+2-g) / ( 3*q+2*(1-g))
                                        #  cat(sprintf("q=%i B=%f A1=%f A2=%f\n", q, B, A1, A2))
                                        #  cat(sprintf("g=%f, c1=%f, c2=%f, c3=%f\n", g, c1, c2, c3))
###orgDef: E<-1/(1-A2/q)
###orgDef: V<- 2/q * (1+c1*B) /  ( (1-c2*B)^2 * (1-c3*B) )
  
  ##EE     <- 1/(1-A2/q)
                                        #VV     <- (2/q) * (1+c1*B) /  ( (1-c2*B)^2 * (1-c3*B) )
  EE     <- 1 + (A2/q)
  VV     <- (2/q)*(1+B)
  EEstar <- 1/(1-A2/q)
  VVstar <- (2/q)*((1+c1*B)/((1-c2*B)^2 * (1-c3*B)))
                                        #  cat(sprintf("EE=%f VV=%f EEstar=%f VVstar=%f\n", EE, VV, EEstar, VVstar))
  
  V0<-1+c1*B
  V1<-1-c2*B
  V2<-1-c3*B
  V0<-ifelse(abs(V0)<1e-10,0,V0)
                                        #  cat(sprintf("V0=%f V1=%f V2=%f\n", V0, V1, V2)) 
  
###orgDef: V<- 2/q* V0 /(V1^2*V2)
###orgDef: rho <-  V/(2*E^2)
  
  rho <- 1/q * (.divZero(1-A2/q,V1))^2 * V0/V2
  df2 <- 4 + (q+2)/ (q*rho-1)
                                        #  cat(sprintf("rho=%f, df2=%f\n", rho, df2))
  
###orgDef: F.scaling <-  df2 /(E*(df2-2))
###altCalc F.scaling<- df2 * .divZero(1-A2/q,df2-2,tol=1e-12)
  ## this does not work because df2-2 can be about 0.1
  F.scaling<-ifelse( abs(df2-2)<1e-2, 1 , df2*(1-A2/q)/(df2-2))
  
### The F-statistic; scaled and unscaled
  betaDiff<-cbind(beta-beta0)
  Wald  <- as.numeric(t(betaDiff) %*% t(L) %*% solve(L%*%PhiA%*%t(L), L%*%betaDiff))
  WaldU <- as.numeric(t(betaDiff) %*% t(L) %*% solve(L%*%Phi%*%t(L), L%*%betaDiff))
                                        #  cat(sprintf("Wald=%f WaldU=%f\n", Wald, WaldU))
  
  Fstat  <- F.scaling/q * Wald
  pval   <- pf(Fstat,df1=q,df2=df2,lower.tail=FALSE)
  FstatU <- Wald/q
  pvalU  <- pf(FstatU,df1=q,df2=df2,lower.tail=FALSE)
  
### SHD addition: calculate bartlett correction and gamma approximation
###
##   ## Bartlett correction - X2 distribution
##   BCval   <- 1 / EE
##   BCstat  <- BCval * Wald
##   p.BC    <- 1-pchisq(BCstat,df=q)
## #  cat(sprintf("Wald=%f BCval=%f BC.stat=%f p.BC=%f\n", Wald, BCval, BCstat, p.BC))
##   ## Gamma distribution
##   scale   <- q*VV/EE
##   shape   <- EE^2/VV
##   p.Ga    <- 1-pgamma(Wald, shape=shape, scale=scale)
## #  cat(sprintf("shape=%f scale=%f p.Ga=%f\n", shape, scale, p.Ga))

  stats<-c(df1=q,df2=df2,Fstat=Fstat,
           p.value=pval, F.scaling=F.scaling, FstatU=FstatU, p.value.U=pvalU,
           A1=A1, A2=A2, V0=V0, V1=V1, V2=V2, rho=rho)
  stats

}


print.KRmodcomp <- function(x,...){

    cat(sprintf("F-test with Kenward-Roger approximation; computing time: %.2f sec.\n",
                attr(x,"ctime")))
    ##     formLarge <- x$f.large
    ##     attributes(formLarge) <- NULL
    cat("Large : ")
    print(x$f.large)
    if (inherits(x$f.small,"call"))
    {
      cat("small : ")
      print(x$f.small)
    }
    else {
      formSmall <- x$f.small
      cat("small : Lbeta=beta0")
      cat('L=')
      print(formSmall$L)
      cat('beta0=')
      print(formSmall$beta0)
    }

    stats<-x$stats
##     cat(sprintf("df1=%3i, df2=%8.2f, Fstat=%8.2f, pval=%7.5f, Fscal= %4.3f \n",
##                 stats['df1'], stats['df2'], stats['Fstat'], stats['pval'],
##                 stats['F.scaling']) )

    sss<-stats[c("Fstat","df1","df2","p.value","F.scaling")]
    sss<-as.data.frame(as.list(sss))
    
    sss$p.value <-round(sss$p.value, options("digits")$digits)
    sss$Fstat   <-round(sss$Fstat,   options("digits")$digits)

    print(sss, row.names=FALSE)
    
##    print(as.data.frame(stats))
    
    if (stats['F.scaling']<0.2) {
        cat('The scaling factor for the F-statistic is smaller than 0.2 \n')
        cat('The unscaled statistic might be more reliable \n ')
        cat('Results fromm the unscaled F-statistic \n')
        cat(sprintf("df1=%3i, df2=%8.2f, FstatU=%8.2f, pvalU=%7.5f  \n",
                    stats['df1'], stats['df2'], stats['FstatU'], stats['pvalU']))
    }
    return(invisible(x))
}
