\name{KRmodcomp_init}
\alias{KRmodcomp_init}
\alias{KRmodcomp_init.mer}
\alias{KRmodcomp_init.lmerMod}

\title{Model input check}
\description{
  Utility functions check the models for KRmodcomp.lmer function.
  Not to be used directly by the user
}
\usage{
KRmodcomp_init(m1, m2, matrixOK=FALSE)
\method{KRmodcomp_init}{mer}(m1, m2, matrixOK=FALSE)
\method{KRmodcomp_init}{lmerMod}(m1, m2, matrixOK=FALSE)
}

\arguments{
  \item{m1}{lmer model}
  \item{m2}{lmer model or restriction matrix}
  \item{matrixOK}{if FALSE then m2 must not be a restriction matrix but a lmer model,
    if TRUE m2 may also be a resteriction matrix}
}
\details{
  The function checks whether m1 is an lmer model
  and m2 either an lmer model or a restriction matrix.
  Furthermore, it is checkes whether the fixed effects space of
  m2 is a subspace of m1.
  In the following <m1> is the fixed efects space of model m1.
}
\value{
  \item{code}{-1: <m1> =<m2>, 0: <m1> < <m2>, 1: <m1> > <m2>}  
}

\author{
  Ulrich Halekoh \email{ulrich.halekoh@agrsci.dk}, Søren Højsgaard \email{sorenh@math.aau.dk}
}

\seealso{
  \code{\link{lmer}}.
}

\keyword{function}
