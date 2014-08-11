\name{get_ddf_Lb}
\alias{get_ddf_Lb}
\alias{get_ddf_Lb.lmerMod}
\alias{ddf_Lb}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Adjusted denomintor degress freedom for linear estimate for linear
 mixed model.
}
\description{
  Get adjusted denomintor degress freedom for testing L'beta=0 in a
  linear mixed model. 
}
\usage{
get_ddf_Lb(object, Lcoef)
ddf_Lb(VVa, Lcoef, VV0=VVa)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A linear mixed model object.
}
  \item{Lcoef}{
    A vector with the same length as \code{fixef(object)}
  }
  \item{VVa, VV0}{
    Adjusted and unadjusted covariance matrix for the fixed effects
    parameters. Adjusted covariance matrix is obtained with
    \code{vcovAdj()} and the unadjusted with \code{vcov()}.
    }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Adjusted degrees of freedom (adjusment made by a Kenward-Roger
  approximation).
}

\references{
  Ulrich Halekoh, Søren Højsgaard (2014).,
  A Kenward-Roger Approximation and Parametric Bootstrap Methods for Tests in Linear Mixed Models - The R Package pbkrtest.,
  Journal of Statistical Software, 58(10), 1-30., \url{http://www.jstatsoft.org/v59/i09/}
}

\author{
  Søren Højsgaard, \email{sorenh@math.aau.dk}
}
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{KRmodcomp}}
  \code{\link{vcovAdj}}
}
\examples{
## To be written
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

