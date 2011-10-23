modcomp_init <- function(m1, m2, matrixOK){
    UseMethod("modcomp_init")
}

modcomp_init.mer<-function(m1,m2,matrixOK=FALSE) {
                                        #comparison of the mean structures of the models
                                        #is is  tested for that
                                        # m1 is mer and
                                        # m2 either mer or a matrix
                                        #checking for mer models
  mers<- if ( 'mer' %in% class(m1) &
             ( ('mer' %in% class(m2)) | is.matrix(m2) ) ) TRUE else FALSE
  
  if (!mers) {
    print('Error in modcomp_init')
    print(paste('either model ',substitute(m1), ' is not of class mer',sep=' '))
    print(paste('or model ', substitute(m2),' is neither of class mer nor matrix',sep=''))
    stop()
  }
  
                                        #checking matrixcOK is FALSE but m2 is a matrix
  if (!matrixOK & is.matrix(m2)) {
    cat ('Error in modcomp_init \n')
    cat (paste('matrixOK =FALSE but the second model: ', substitute(m2),
               '\n is  specified via a restriction matrix \n \n',sep=''))
    stop()
  }
  
  
  Xlarge<-m1@X
  rlarge<-rankMatrix(Xlarge)
  code<- if ('mer' %in% class(m2)) {
    Xsmall<-m2@X
    rsmall<-rankMatrix(Xsmall)
    rboth<- rankMatrix(cbind(Xlarge,Xsmall))
    if (rboth == pmax(rlarge,rsmall)) {
      if (rsmall< rlarge) {1} else {
        if (rsmall > rlarge) 0 else -1
      }
    }
    else {-1}
  } else {
                                        #now model m2  is a restriction matrix
    if (rankMatrix(rbind(Xlarge,m2)) > rlarge) -1 else 1
  }
  code
}


