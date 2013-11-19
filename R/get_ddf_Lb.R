get_ddf_Lb <- function(object, Lcoef){
    UseMethod("get_ddf_Lb")
}

get_ddf_Lb.lmerMod <- function(object, Lcoef){
    ddf_Lb(vcovAdj(object), Lcoef)
}

ddf_Lb <- function(VVa, Lcoef, VV0=VVa){

    .spur = function(U){
        sum(diag(U))
    }
    .divZero = function(x,y,tol=1e-14){
        ## ratio x/y is set to 1 if both |x| and |y| are below tol
        x.y  =  if( abs(x)<tol & abs(y)<tol) {1} else x/y
        x.y
    }

    vlb = sum(Lcoef * (VV0 %*% Lcoef))
    Theta = Matrix(as.numeric(outer(Lcoef, Lcoef) / vlb), nrow=length(Lcoef))
    P = attr(VVa, "P")
    W = attr(VVa, "W")

    A1 = A2 = 0
    ThetaVV0 = Theta%*%VV0
    n.ggamma = length(P)
    for (ii in 1:n.ggamma) {
        for (jj in c(ii:n.ggamma)) {
            e = ifelse(ii==jj, 1, 2)
            ui = ThetaVV0 %*% P[[ii]] %*% VV0
            uj = ThetaVV0 %*% P[[jj]] %*% VV0
            A1 =  A1 +  e* W[ii,jj] * (.spur(ui) * .spur(uj))
            A2 =  A2 +  e* W[ii,jj] *  sum(ui * t(uj))
        }}

    ## substituted q = 1 in pbkrtest code and simplified
    B  =  (A1 + 6 * A2) / 2
    g  =  (2 * A1 - 5 * A2)  / (3 * A2)
    c1 =  g/(3 + 2 * (1 - g))
    c2 =  (1 - g) / (3 + 2 * (1 - g))
    c3 =  (3 - g) / (3 + 2 * (1 - g))
    EE =  1 + A2
    VV =  2 * (1 + B)
    EEstar  =  1/(1 - A2)
    VVstar  =  2 * ((1 + c1 * B)/((1 - c2 * B)^2  *  (1 - c3 * B)))
    V0 = 1 + c1 * B
    V1 = 1 - c2 * B
    V2 = 1 - c3 * B
    V0 = ifelse(abs(V0) < 1e-10, 0, V0)
    rho  = (.divZero(1 - A2, V1))^2 * V0/V2
    df2  =  4 + 3 / (rho - 1)
    ## cat(sprintf("Lcoef: %s\n", toString(Lcoef)))
    ## cat(sprintf("df2: %f\n", df2))
    df2
}

## .get_ddf: Adapted from Russ Lenths 'lsmeans' package.
## Returns denom d.f. for testing lcoefs'beta = 0 where lcoefs is a vector
# VVA is result of call to VVA = vcovAdj(object, 0) in pbkrtest package
# VV is vcov(object) ## May not now be needed
# lcoefs is contrast of interest
# varlb is my already-computed value of lcoef' VV lcoef = est variance of lcoef'betahat


.get_ddf <- function(VVa, VV0, Lcoef, varlb) {

    ## print(VVa); print(VV0)
    ##   ss<<-list(VVa=VVa, VV0=VV0)

  .spur = function(U){
      ##print(U)
    sum(diag(U))
  }
  .divZero = function(x,y,tol=1e-14){
    ## ratio x/y is set to 1 if both |x| and |y| are below tol
    x.y  =  if( abs(x)<tol & abs(y)<tol) {1} else x/y
    x.y
  }

  vlb = sum(Lcoef * (VV0 %*% Lcoef))
  Theta = Matrix(as.numeric(outer(Lcoef, Lcoef) / vlb), nrow=length(Lcoef))
  P = attr(VVa, "P")
  W = attr(VVa, "W")

  A1 = A2 = 0
  ThetaVV0 = Theta%*%VV0
  n.ggamma = length(P)
  for (ii in 1:n.ggamma) {
    for (jj in c(ii:n.ggamma)) {
      e = ifelse(ii==jj, 1, 2)
      ui = ThetaVV0 %*% P[[ii]] %*% VV0
      uj = ThetaVV0 %*% P[[jj]] %*% VV0
      A1 =  A1 +  e* W[ii,jj] * (.spur(ui) * .spur(uj))
      A2 =  A2 +  e* W[ii,jj] *  sum(ui * t(uj))
    }}

  ## substituted q = 1 in pbkrtest code and simplified
  B  =  (A1 + 6 * A2) / 2
  g  =  (2 * A1 - 5 * A2)  / (3 * A2)
  c1 =  g/(3 + 2 * (1 - g))
  c2 =  (1 - g) / (3 + 2 * (1 - g))
  c3 =  (3 - g) / (3 + 2 * (1 - g))
  EE =  1 + A2
  VV =  2 * (1 + B)
  EEstar  =  1/(1 - A2)
  VVstar  =  2 * ((1 + c1 * B)/((1 - c2 * B)^2  *  (1 - c3 * B)))
  V0 = 1 + c1 * B
  V1 = 1 - c2 * B
  V2 = 1 - c3 * B
  V0 = ifelse(abs(V0) < 1e-10, 0, V0)
  rho  = (.divZero(1 - A2, V1))^2 * V0/V2
  df2  =  4 + 3 / (rho - 1)
  ## cat(sprintf("Lcoef: %s\n", toString(Lcoef)))
  ## cat(sprintf("df2: %f\n", df2))

  df2
}
