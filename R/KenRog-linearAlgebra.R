.spur<-function(U){
  sum(diag(U))
}

.spurAB<-function(A,B){
  sum(A*t.default(B))
}
# if A eller B er symmetrisk så er trace(AB)=sum(A*B)

.fatAB<-function(A,B) {
  ## <A> in <B>
  ## determine L such that  <A>={Bb| b in Lb=0}
  
  d<-qr(cbind(A,B))$rank - qr(B)$rank
  if (d>0) {
    print('Error:  <A> not subspace of <B> ')
    return()
  }
  k<-qr(cbind(A,B))
  k<-qr.Q(k)[,1:k$rank]
  L<-(ginv(k) %*% B)[-(1:qr(A)$rank),,drop=FALSE]
                                        # making th ros of L orthogonal
  L<-t(qr.Q(qr(t(L))))
  L<-ifelse(abs(L)<1e-15,0,L)
  L
}

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
  rW<-rankMatrix(W)
  k<-qr(cbind(W,diag(nrow(W))  ) )
  Worth<-qr.Q(k)[,c( (rW+1):k$rank),drop=FALSE]
  Worth
}

