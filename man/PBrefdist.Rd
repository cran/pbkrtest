\name{PBrefdist}
\alias{PBrefdist}
\alias{PBrefdist.mer}
\alias{PBrefdist.merMod}
\alias{PBrefdist.lm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate reference distribution using parametric bootstrap
}
\description{
    Calculate reference distribution of likelihood ratio statistic in
    mixed effects models using parametric bootstrap 
}
\usage{
PBrefdist(largeModel, smallModel, nsim = 1000, seed=NULL, cl = NULL, details = 0)
%\method{PBrefdist}{lm}(largeModel, smallModel, nsim = 1000, seed=NULL, cl = NULL, details = 0)
%\method{PBrefdist}{mer}(largeModel, smallModel, nsim = 1000, seed=NULL, cl = NULL, details = 0)
%\method{PBrefdist}{merMod}(largeModel, smallModel, nsim = 1000, seed=NULL, cl = NULL, details = 0)
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
  \item{seed}{
    Seed for the random number generation.
  }
  \item{cl}{
    A vector identifying a cluster; used for
    calculating the reference distribution using several cores. See
    examples below.
  }
  \item{details}{
    The amount of output produced. Mainly relevant for debugging
    purposes.
  }
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
A numeric vector
}
%\references{
%% ~put references to the literature/web site here ~
%}

\references{
  Ulrich Halekoh, Søren Højsgaard (2014).,
  A Kenward-Roger Approximation and Parametric Bootstrap Methods for Tests in Linear Mixed Models - The R Package pbkrtest.,
  Journal of Statistical Software, 58(10), 1-30., \url{http://www.jstatsoft.org/v59/i09/}
}


\author{
  Soren Hojsgaard \email{sorenh@math.aau.dk}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{PBmodcomp}},
%  \code{\link{BCmodcomp}},
%  \code{\link{makePBcluster}},
%  \code{\link{stopPBcluster}},
  \code{\link{KRmodcomp}}
}
\examples{
data(beets)
head(beets)
beet0<-lmer(sugpct~block+sow+harvest+(1|block:harvest), data=beets, REML=FALSE)
beet_no.harv <- update(beet0, .~.-harvest)
rr <- PBrefdist(beet0, beet_no.harv, nsim=20)
rr

## Note clearly many more than 10 simulations must be made in practice.

## Computations can be made in parallel using several processors:
\dontrun{
cl <- makeSOCKcluster(rep("localhost", 4))
clusterEvalQ(cl, library(lme4))
clusterSetupSPRNG(cl)
rr <- PBrefdist(beet0, beet_no.harv, nsim=20)
stopCluster(cl)
}
## Above, 4 cpu's are used and 5 simulations are made on each cpu.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{models}
