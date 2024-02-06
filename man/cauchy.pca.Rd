\name{Robust PCA using the Cauchy distribution}
\alias{cauchy.pca}
\title{
Robust PCA using the Cauchy distribution
}

\description{
Robust PCA using the Cauchy distribution.
}

\usage{
cauchy.pca(x, k = 1, center = "sm", scale = "mad", trials = 20, parallel = FALSE)
}

\arguments{
\item{x}{
A numerical matrix with the data.
}
\item{k}{
The number of eigenvectors to extract.
}
\item{center}{
The way to center the data. This can be either "sm" corresponding to the spatial median,
"med" corresponding to the classical variable-wise median. Alternatively the user can specify their own vector.
}
\item{scale}{
This is the method to scale the data. The default value is "mad" corresponding to the
mean absolute deviation, computed column-wise. Alternatively the user can provide their own
vector.
}
\item{trials}{
The number of trials to attempt. How many times the algorithm will be performed with different starting values (different starting vectors).
}
\item{parallel}{
If you want parallel computations set this equal to TRUE.
}
}

\details{
This is the main function used to extract the Cauchy robust eigenvectors.
}

\value{
A list including:
\item{runtime}{
The duration (in seconds) of the algorithm.
}
\item{loglik}{
The minimum maximum Cauchy log-likelihood.
}
\item{mu}{
The estimated location parameter of the Cauchy ditribution.
}
\item{su}{
The estimated scale parameter of the Cauchy ditribution.
}
\item{loadings}{
A matrix with the robust eigenvectors.
}
}

\references{
Fayomi A., Pantazis Y., Tsagris M. and Wood A.T.A. (2024). Cauchy robust principal component analysis with applications to high-dimensional data sets. Statistics and Computing, 34: 26.
https://doi.org/10.1007/s11222-023-10328-x
}

\author{
Michail Tsagris, Aisha Fayomi, Yannis Pantazis and Andrew T.A. Wood.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}


\seealso{
\code{ \link{cauchy.mle}
}
}

\examples{
x <- as.matrix( iris[, 1:4] )
cauchy.pca(x, k = 1)
}



