\name{vcovAdj}
\alias{vcovAdj}
\title{Ajusted covariance matrix for linear mixed models according to Kenward and Roger}
\description{
  Kenward and Roger (1997) describbe an improved small sample approximation to
  the covariance matrix estimate of the fixed parameters in a linear mixed model.
}
\usage{
vcovAdj(largeModel, details=0)
}
	
\arguments{
  \item{largeModel}{An \code{lmer} model}
  \item{details}{If larger than 0 some timing details are printed.}
}

\value{
  \item{phiA}{the estimated covariance matrix, this has attributed P, a list of matrices used in \code{KR_adjust}
    and the estimated matrix  W of the variances of the covariance parameters of the random effetcs}
}

\references{
  Kenward, M. G. and Roger, J. H. (1997), \emph{Small Sample Inference 
    for Fixed Effects from Restricted Maximum Likelihood}, Biometrics 53: 983-997. }

\note{
  This functionality is not thoroughly tested and should be used with
  care. Please do report bugs etc.
}

\author{
  Ulrich Halekoh \email{ulrich.halekoh@agrsci.dk}, Søren Højsgaard \email{sorenh@math.aau.dk}
}


\seealso{
  \code{\link{lmer}}.
}
\examples{
(fmLarge <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
## removing Day
(vcovAdj(fmLarge,detail=0))
}
\keyword{function}
