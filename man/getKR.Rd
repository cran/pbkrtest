\name{getKR}
\alias{getKR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extract (or "get") components from a \code{KRmodcomp} object.
}
\description{
  Extract (or "get") components from a \code{KRmodcomp} object, which is
  the result of the \code{KRmodcomp} function. 
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getKR(object, name = c("ndf", "ddf", "Fstat", "p.value",
     "F.scaling", "FstatU", "p.valueU", "aux"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A \code{KRmodcomp} object, which is
    the result of the \code{KRmodcomp} function
  }
  \item{name}{
    The available slots. If \code{name} is missing or \code{NULL} then
    everything is returned. 
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
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
  \code{\link{KRmodcomp}}
  \code{\link{PBmodcomp}}  
  \code{\link{vcovAdj}}  
}
\examples{
data(beets, package='pbkrtest')
lg <- lmer(sugpct ~ block + sow + harvest + (1|block:harvest), 
              data=beets, REML=FALSE)
sm <- update(lg, .~. - harvest)
xx<-KRmodcomp(lg, sm)
getKR(xx, "ddf") # get denominator degrees of freedom.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

