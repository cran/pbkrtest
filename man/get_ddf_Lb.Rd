\name{get_Lb_ddf}
\alias{get_Lb_ddf}
\alias{get_Lb_ddf.lmerMod}
\alias{Lb_ddf}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Adjusted denomintor degress freedom for linear estimate for linear
 mixed model.
}
\description{
  Get adjusted denomintor degress freedom for testing Lb=0 in a
  linear mixed model where L is a restriction matrix.
}
\usage{
get_Lb_ddf(object, L)
Lb_ddf(L, V0, Vadj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A linear mixed model object.
  }

  \item{L}{
    A vector with the same length as \code{fixef(object)} or a matrix
    with the same number of columns as the length of \code{fixef(object)}
  }

  \item{V0, Vadj}{
    Unadjusted and adjusted covariance matrix for the fixed effects
    parameters. Undjusted covariance matrix is obtained with
    \code{vcov()} and adjusted with \code{vcovAdj()}.
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
  Soren Hojsgaard, \email{sorenh@math.aau.dk}
}
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{KRmodcomp}},
  \code{\link{vcovAdj}},
  \code{\link{model2restrictionMatrix}},
  \code{\link{restrictionMatrix2model}}

}
\examples{
(fmLarge <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
## removing Days
(fmSmall <- lmer(Reaction ~ 1 + (Days|Subject), sleepstudy))
anova(fmLarge,fmSmall)

KRmodcomp(fmLarge,fmSmall)  ## 17 denominator df's
get_Lb_ddf(fmLarge, c(0,1)) ## 17 denominator df's

# Notice: The restriction matrix L corresponding to the test above
# can be found with
L<-model2restrictionMatrix(fmLarge, fmSmall)
L
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

