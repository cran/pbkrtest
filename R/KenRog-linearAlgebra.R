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
# let QR be the QR desoposition of matrix A|B
# whre A|B means horizontal concatenation
# such that Q=(Q1,Q2) where
# <Q1>= <A>
# then Bb= Q1 x1 + Q2 x2
# for some x1 and x2 and
#  <A> = {Q1 x1 + Q2 x2 | x2=0}
# since
#Q1=    
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

