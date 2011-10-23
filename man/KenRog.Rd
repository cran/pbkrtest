\name{KenwardRoger}
\alias{KRmodcomp}
\alias{KRmodcomp.mer}
\title{Ftest and degrees of freedom based on Kenward-Roger approximation}
\description{
  An approximate F-test based on the Kenward-Roger 
  approach.
}
\usage{
KRmodcomp(largeModel, smallModel, beta0=0)
\method{KRmodcomp}{mer}(largeModel, smallModel, beta0=0)
}
	
\arguments{
  \item{largeModel}{A  lmer model}
  \item{smallModel}{A  lmer model or a restriction matrix}
  \item{beta0}{A number or a vector of the beta of the hypothesis, 
    e.g. L beta=beta0. beta0=0 if modelSmall is a model not a
    restriction matrix.}
%  \item{x}{A KRmodcomp object}
  \item{\dots}{Additonal arguments to print function}	
}

\details{
  An F test is calculated according to the approach of Kenward
  and Roger (1997).  The function works for linear mixed models fitted
  with the \code{lmer} function of the \pkg{lme4} package. Only models
  where the covariance structure is a sum of known matrices can be
  compared.

  The \code{largeModel} may be a model fitted with \code{lmer}
  either using \code{REML=TRUE} or \code{REML=FALSE}.  The
  \code{smallModel} can be a model fitted with \code{lmer}. It must have
  the same coavarianc structure as \code{largeModel}. Furthermore, its
  linear space of expectation must be a subspace of the space for
  \code{largeModel}.  The model \code{smallModel} can also be a
  restriction matrix \code{L} specifying the hypothesis \code{L beta
  =beta0}, where \code{L} is a \eqn{k \times p}{k X p} matrix and
  \eqn{\beta} is a \eqn{p} column vector the same length as
  \code{fixef(largeModel)}.

  \eqn{\beta0} is a \eqn{k} column vector.

  Notice: It cannot be guaranteed that the results agree with other
  implementations of the Kenward-Roger approach!
}

\references{
  Kenward, M. G. and Roger, J. H. (1997), \emph{Small Sample Inference 
    for Fixed Effects from Restricted Maximum Likelihood}, Biometrics 
  53: 983-997. }

\note{
  This functionality is not thoroughly tested and should be used with
  care. Please do report bugs etc.
}

\seealso{
  \code{\link{lmer}}.
}
\examples{
(fmLarge <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
## removing Day
(fmSmall <- lmer(Reaction ~ 1 + (Days|Subject), sleepstudy))
(KRmodcomp(fmLarge,fmSmall))

## The same test using a restriction matrix
L<-cbind(0,1)
(KRmodcomp(fmLarge,L))
}
\keyword{function}
