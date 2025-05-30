% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SAT_modcomp.R
\name{get_Fstat_ddf}
\alias{get_Fstat_ddf}
\title{Compute denominator degrees of freedom for F-test}
\usage{
get_Fstat_ddf(nu, tol = 1e-08)
}
\arguments{
\item{nu}{vector of denominator degrees of freedom for the
t-statistics}

\item{tol}{tolerance on the consecutive differences between
elements of nu to determine if mean(nu) should be returned}
}
\value{
the denominator degrees of freedom; a numerical scalar
}
\description{
From a vector of denominator degrees of freedom from independent
t-statistics (\code{nu}), the denominator degrees of freedom for
the corresponding F-test is computed.
}
\details{
Note that if any \code{nu <= 2} then \code{2} is returned. Also, if
all nu are within \code{tol} of each other the simple average of the
nu-vector is returned.  This is to avoid downward bias.
}
\author{
Rune Haubo B. Christensen. Adapted to pbkrtest by Søren Højsgaard.
}
\keyword{internal}
