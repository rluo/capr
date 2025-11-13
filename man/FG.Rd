\name{FG}
\alias{FG}
\title{Flury-Gautschi Common Principal Components}
\usage{
FG(cov_array, P = NULL, M = NULL, maxit = 30L)
}
\arguments{
\item{cov_array}{Numeric array of dimension \eqn{P \times P \times M}
containing covariance matrices in its slices.}

\item{P}{Optional integer specifying the dimension of each covariance
matrix. Defaults to \code{dim(cov_array)[1]}.}

\item{M}{Optional integer specifying the number of matrices/slices.
Defaults to \code{dim(cov_array)[3]}.}

\item{maxit}{Integer scalar; number of outer iterations.}
}
\value{
A \eqn{P \times P} numeric matrix containing the estimated common
loadings.
}
\description{
Implements the Flury and Gautschi (1986) algorithm to iteratively align a
collection of covariance matrices stored in a 3D array to a shared loading
structure by sweeping over all pairs of variable indices and updating local
rotations.
}
\examples{
\dontrun{
set.seed(1)
P <- 3
M <- 4
mats <- replicate(M, {
    A <- matrix(rnorm(P * P), P, P)
    crossprod(A)
}, simplify = FALSE)
cov_cube <- array(NA_real_, dim = c(P, P, M))
for (k in 1:M) cov_cube[, , k] <- mats[[k]]
FG(cov_cube, maxit = 5)
}
}
