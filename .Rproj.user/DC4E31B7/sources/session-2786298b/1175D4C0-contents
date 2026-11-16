\name{MLE of the Cauchy distribution}
\alias{cauchy.mle}
\title{
MLE of the Cauchy distribution
}
\description{
MLE of the Cauchy distribution.
}

\usage{
cauchy.mle(x, tol = 1e-07)
}

\arguments{
\item{x}{
A numerical vector with data.
}
\item{tol}{
The tolerance level up to which the maximisation stops set to 1e-09 by default.
}
}

\details{
Instead of maximising the log-likelihood via a numerical optimiser we have used a Newton-Raphson algorithm which 
is faster. The Cauchy is the t distribution with 1 degree of freedom. 
}

\value{
A list including: 
\item{iters}{The number of iterations required for the Newton-Raphson to converge.}
\item{loglik}{The value of the maximised log-likelihood.}
\item{param}{The vector of the parameters.}
}

\references{
Johnson, Norman L. Kemp, Adrianne W. Kotz, Samuel (2005). Univariate Discrete Distributions (third edition). 
Hoboken, NJ: Wiley-Interscience.

https://en.wikipedia.org/wiki/Wigner_semicircle_distribution
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{cauchy.pca}
}
}

\examples{
x <- rcauchy(1000)
a <- cauchy.mle(x)
}
