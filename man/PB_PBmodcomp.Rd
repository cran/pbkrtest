\name{PBmodcomp}
\alias{PBmodcomp}
\alias{PBmodcomp.mer}
%\alias{BCmodcomp}
%\alias{BCmodcomp.mer}
\alias{makePBcluster}
\alias{stopPBcluster}
\alias{plot.XXmodcomp}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Model comparison of mixed models using parametric bootstrap methods.
}
\description{
  Model comparison of mixed models using parametric bootstrap methods.
}
\usage{
%BCmodcomp(largeModel, smallModel, nsim = 200, ref = NULL, cl = NULL, details = 0)
%\method{BCmodcomp}{mer}(largeModel, smallModel, nsim = 200, ref = NULL, cl = NULL, details = 0)
PBmodcomp(largeModel, smallModel, nsim = 200, ref = NULL, cl = NULL, details = 0)
\method{PBmodcomp}{mer}(largeModel, smallModel, nsim = 200, ref = NULL, cl = NULL, details = 0)
makePBcluster(n=2)
stopPBcluster()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{largeModel}{
    A linear mixed effects model as fitted with the \code{lmer()}
    function in the \pkg{lme4} package. This model muse be larger than
    \code{smallModel} (see below).
  }
  \item{smallModel}{
    A linear mixed effects model as fitted with the \code{lmer()}
    function in the \pkg{lme4} package. This model muse be smaller than
    \code{largeModel} (see above).
  }
  \item{nsim}{
    The number of simulations to form the reference distribution.
  }
  \item{ref}{
    Vector containing samples from the reference distribution. If NULL,
    this vector will be generated using PBrefdist(). }
  \item{cl}{
    A vector identifying a cluster; used for
    calculating the reference distribution using several cores. See
    examples below.
  }
  \item{n}{
    Number of clusters to be used in computations; a typical choice of
    \code{n} would be the number of CPUs available.
  }
  \item{details}{
    The amount of output produced. Mainly relevant for debugging
    purposes.
  }
}
\details{
  Under the fitted hypothesis (i.e. under the fitted small model)
  \code{nsim} samples of the likelihood ratio test statistic (LRT) are
  generetated.

  Then p-values are calculated as follows:

  LRT: Assuming that LRT has a chi-square distribution.

  PBtest: The fraction of simulated LRT-values that are larger or equal
  to the observed LRT value.

  Bartlett: A Bartlett correction is of LRT is calculated from the mean
  of the simulated LRT-values

  Gamma: The reference distribution of LRT is assumed to be a gamma
  distribution with mean and variance determined as the sample mean and
  sample variance of the simulated LRT-values.

}
%\value{
%  ## To be added...
%}
%\references{
%  ## To be added...
%}
\author{
  Søren Højsgaard <sorenh at mail dot dk>
}
\note{
  This functionality is not thoroughly tested and should be used with
  care. Please do report bugs etc.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{KRmodcomp}}
  \code{\link{PBrefdist}}
}
\examples{
data(beets)
head(beets)
beet0<-lmer(sugpct~block+sow+harvest+(1|block:harvest), data=beets, REML=FALSE)
beet_no.harv <- update(beet0, .~.-harvest)
PBmodcomp(beet0, beet_no.harv, nsim=20)

\dontrun{
## Vanilla
PBmodcomp(beet0, beet_no.harv)
## Simulate reference distribution separately:
rr <- PBrefdist(beet0, beet_no.harv, nsim=20)
PBmodcomp(beet0, beet_no.harv, ref=rr)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{models}% __ONLY ONE__ keyword per line
\keyword{inference}