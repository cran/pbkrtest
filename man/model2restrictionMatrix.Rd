\name{model2restrictionMatrix}
\alias{model2restrictionMatrix}
\alias{model2restrictionMatrix.lm}
\alias{model2restrictionMatrix.mer}
\alias{model2restrictionMatrix.merMod}
\alias{restrictionMatrix2model}
\alias{restrictionMatrix2model.lm}
\alias{restrictionMatrix2model.mer}
\alias{restrictionMatrix2model.merMod}

\title{
  Conversion between a model object and a restriction matrix
}
\description{
  Testing a small model under a large model corresponds imposing
  restrictions on the model matrix of the larger model and these
  restrictions come in the form of a restriction matrix. These functions
  converts a model to a restriction matrix and vice versa.
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
model2restrictionMatrix(largeModel, smallModel)
restrictionMatrix2model(largeModel, LL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{largeModel, smallModel}{
    Model objects of the same "type". Possible types are linear mixed
    effects models and linear models (including generalized linear models)
  }
  \item{LL}{
    A restriction matrix.
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  model2restrictionMatrix: A restriction matrix.

  restrictionMatrix2model: A model object.
}

\references{
  Ulrich Halekoh, Søren Højsgaard (2014).,
  A Kenward-Roger Approximation and Parametric Bootstrap Methods for Tests in Linear Mixed Models - The R Package pbkrtest.,
  Journal of Statistical Software, 58(10), 1-30., \url{http://www.jstatsoft.org/v59/i09/}
}


\author{
  Ulrich Halekoh \email{uhalekoh@health.sdu.dk}, Soren Hojsgaard \email{sorenh@math.aau.dk}
}

\note{
  That these functions are visible is a recent addition; minor changes
  may occur.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{PBmodcomp}},  \code{\link{PBrefdist}},  \code{\link{KRmodcomp}}
}
\examples{
library(pbkrtest)
data("beets", package = "pbkrtest")
sug <- lm(sugpct ~ block + sow + harvest, data=beets)
sug.h <- update(sug, .~. - harvest)
sug.s <- update(sug, .~. - sow)

## Construct restriction matrices from models
L.h <- model2restrictionMatrix(sug, sug.h); L.h
L.s <- model2restrictionMatrix(sug, sug.s); L.s

## Construct submodels from restriction matrices
mod.h <- restrictionMatrix2model(sug, L.h); mod.h
mod.s <- restrictionMatrix2model(sug, L.s); mod.s

## The models have the same fitted values
plot(fitted(mod.h), fitted(sug.h))
plot(fitted(mod.s), fitted(sug.s))
## and the same log likelihood
logLik(mod.h)
logLik(sug.h)
logLik(mod.s)
logLik(sug.s)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{models}% __ONLY ONE__ keyword per line
