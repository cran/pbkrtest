.spur<-function(U){
  sum(diag(U))
}

.spurAB<-function(A,B){
  sum(A*t.default(B))
}
# if A eller B er symmetrisk så er trace(AB)=sum(A*B)

.fatBL<-function(B,L) {
  ## find A such that
  ## <A>={Bb| b in Lb=0}
  if ( ncol(B) != ncol(L) ) {
    cat('Error \n number of columns of B and L unequal \n')
    return()
  }
  A<-B %*% .orthComplement(t(L))
  A
}

.orthComplement<-function(W) {
  ##orthogonal complement of <W>: <W>orth=<Worth>
  rW <-rankMatrix(W)
  k  <-qr(cbind(W,diag(nrow(W))  ) )
  Worth<-qr.Q(k)[,c( (rW+1):k$rank),drop=FALSE]
  Worth
}


.restrictionMatrixBA<-function(B,A) {
  ## <A> in <B>
  ## determine L such that  <A>={Bb| b in Lb=0}
  d<-rankMatrix(cbind(A,B)) - rankMatrix(B)
  if (d>0) {
    print('Error:  <A> not subspace of <B> ')
    return()
  }
  Q <-qr.Q(qr(cbind(A,B)))
  Q2<-Q[,(rankMatrix(A)+1):rankMatrix(B)] 
  L <- t(Q2)  %*% B
  ##make rows of L2 orthogonal
  L<-t(qr.Q(qr(t(L))))
L
}


.createRestrictionMatrix <- function (largeModel, smallModel) {
  L<- if(is.matrix(smallModel)) {
    ## ensures  that L is of full row rank:
    q<-rankMatrix(L)
    if (q < nrow(L) ){
      t(qr.Q(qr(t(L)))[,1:qr(L)$rank])
    } else {
      smallModel
    }
  } else  { #smallModel is mer model
    .restrictionMatrixBA(getME(largeModel,'X'),getME(smallModel,'X'))
  }
  L<-.makeSparse(L)
  L
}
