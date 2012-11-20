LMM_Sigma_G  <- function(object, details=0) {
  ##output VAR(Y) =Sigma and the G  matrices

  if (!.is.lmm(object)){
    stop("'object' is not Gaussian linear mixed model")
  }
  
  DB <- details>0
  GGamma <- VarCorr(object)  
  Nindex <- .get_indices(object)
  
  ## number of random effects in each groupFac
  ## residual error here excluded!
  n.groupFac<-Nindex$n.groupFac
  ##  the number of random effects for each grouping factor
  nn.groupFacLevels<-Nindex$nn.groupFacLevels
  
  ## size of the symmetric variance Gamma_i for reach groupFac
  nn.GGamma <- Nindex$nn.GGamma
  ## number of variance parameters of each GGamma_i
  mm.GGamma <-  Nindex$mm.GGamma
  group.index <- Nindex$group.index
 
  ##writing the covariance parameters for the rasndom efefcts in a vector: 
  ggamma <- NULL
  for (ii in 1: (n.groupFac)) {
    Lii<-GGamma[[ii]]
    nu<-ncol(Lii)
    ## The lower.tri construxtion esnures, that (because Lii is symmetric!)
    ## Lii[lower.tri(Lii,diag=TRUE)= Lii[1,1],Lii[1,2],Lii[1,3]..Lii[1,nu],
    ##                               Lii[2,2], Lii[2,3] ...
    ggamma<-c(ggamma,Lii[lower.tri(Lii,diag=TRUE)])
  }
  

  ##adding the residuals variance to ggamma
  ##so in ggamma nd n.ggamma the residual variance is included!
  ggamma<-c(ggamma,sigma(object)^2)
  n.ggamma<-length(ggamma)
  
  Zt<-getME(object,"Zt")
  ## G_r:
  
  t0 <- proc.time()
  
  G<-NULL
  for (ss in 1:n.groupFac)
    {
	ZZ<- .get_Ztgroup(ss, Zt, object) 
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
  nobs <-nrow(getME(object,'X'))
  G    <-c(G,list(sparseMatrix(1:nobs,1:nobs,x=1))) ## The last one is for the residual!
  if(DB){cat(sprintf("Finding G  %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}
  
  Sigma<-ggamma[1]*G[[1]]
  for (ii in 2:n.ggamma) {
    Sigma<- Sigma + ggamma[ii] * G[[ii]]
  }
  
  if(DB){cat(sprintf("Finding Sigma:    %10.5f\n", (proc.time()-t0)[1] )); t0 <- proc.time()}
  
  SigmaG<-list(Sigma=Sigma, G=G, n.ggamma=n.ggamma)
  SigmaG
}  


##### calculation ov adjustes covariance matrix
.get_Nrandomeffects_group <- function(object) {
  ##  \item{object}{A mer or lmerMod  model}
  ## output: number of random effects for each grouping factor
  .cc <- class(object)
  qq<- if (.cc %in% "mer") {
    sapply(object@ST,function(X) nrow(X)) 
  } else {
    sapply(object@cnms, length)
  }
  qq
}
	  
.get_indices <-function(object) {
  ##  \item{object}{A mer or lmerMod  model}
  ## ff= number of random effects terms (..|..)
  ## (..|F1) + (..|F1) are group factors!
  ## without the residual variance
  ##output: list of several inidces

  ff<- getME(object,"n_rtrms")
  
  gg <- sapply(getME(object,"flist"),
             function(x)length(levels(x)))
  
  qq<-.get_Nrandomeffects_group(object)
  
  ## number of variance parameters of each GGamma_i
  ss <- qq * (qq+1)/2

  group.index<-getME(object,"Gp")
			
  list(n.groupFac=ff,
       nn.groupFacLevels=gg,
       nn.GGamma=qq,
       mm.GGamma=ss,
       group.index=group.index)
}

.get_Ztgroup <- function(ii.group,Zt,object) {
  ##arguiments
  ## \item{ii.group}{the index number of a grouping factor}
  ##  \item{Zt}{the transpose of the random factors design matrix Z}
  ##  \item{object}{A mer or lmerMod  model}
  ##output
  ##  submatrix of Zt belongig to grouping factor ii.group
  Nindex <- .get_indices(object)
  nn.groupFacLevels<-Nindex$nn.groupFacLevels
  nn.GGamma <- Nindex$nn.GGamma
  group.index <- Nindex$group.index
  .cc <- class(object)

  zIndex.sub<- if (.cc %in% "mer") {
    Nindex$group.index[ii.group]+
      1+c(0:(nn.GGamma[ii.group]-1))*nn.groupFacLevels[ii.group] +
        rep(0:(nn.groupFacLevels[ii.group]-1),each=nn.GGamma[ii.group]) 
  } else {
    if (.cc %in% "lmerMod" ) {
      c((group.index[ii.group]+1) : group.index[ii.group+1])
    }
  }
  ZZ<-Zt[zIndex.sub,]
  return(ZZ)
}
