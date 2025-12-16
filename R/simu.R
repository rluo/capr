#' Simulate Covariate-Assisted Principal Regression Data
#'
#' Generates synthetic data for the CAP (Covariate Assisted Principal Regression)
#' model. Creates a collection of covariance matrices whose log-spectra vary
#' linearly with binary covariates.
#'
#' @param seed Integer scalar, random seed for reproducibility (default 123).
#' @param p Integer scalar, dimension of the covariance matrices (default 5).
#' @param n Integer scalar, number of samples/observations (default 120).
#'
#' @return A list containing:
#'   \item{S}{Numeric 3D array of dimension \eqn{p \times p \times n};
#'            stack of covariance matrices.}
#'   \item{X}{Numeric matrix of dimension \eqn{n \times 2}; design matrix
#'            with intercept column and binary covariate.}
#'   \item{Q}{Numeric matrix of dimension \eqn{p \times p};
#'            orthogonal matrix (eigenvectors of reference covariance).}
#'   \item{BetaMat}{Numeric matrix of dimension \eqn{2 \times p};
#'                  true regression coefficients of log-eigenvalues.}
#'   \item{H}{Numeric matrix of dimension \eqn{p \times p};
#'            mean of covariance matrices across all samples.}
#'   \item{p}{Integer; dimension parameter.}
#'   \item{n}{Integer; sample size parameter.}
#'
#' @details
#' The generative model:
#' * \eqn{\lambda_i = \exp(X_i \beta)} where \eqn{\beta} are the true log-eigenvalues.
#' * \eqn{S_i = Q \Lambda_i Q^T}, where \eqn{\Lambda_i = \text{diag}(\lambda_i)}.
#' * The true \code{BetaMat} is a fixed \eqn{2 \times p} matrix.
#' * Covariates \eqn{X} are: intercept + Bernoulli(0.5).
#'
#' @examples
#' dat <- simu.capr(seed = 42, p = 5, n = 100)
#' str(dat)
#' # S is a 5 x 5 x 100 array of covariance matrices
#' # X is a 100 x 2 design matrix
#'
#' @export
simu.capr <- function(seed = 123L, p = 5L, n = 120L) {
    stopifnot(length(seed) == 1L, length(p) == 1L, length(n) == 1L)
    set.seed(as.integer(seed))

    BetaMat <- rbind(
        c(5, 4, 1, -1, -2),
        c(0, -1, 1, 0, 0)
    )
    if (ncol(BetaMat) != p) {
        stop("BetaMat must have ", p, " columns.", call. = FALSE)
    }

    Q <- qr.Q(qr(matrix(stats::rnorm(p * p), p, p)))
    X <- cbind(1, stats::rbinom(n, size = 1L, prob = 0.5))

    S <- array(0, dim = c(p, p, n))
    for (i in seq_len(n)) {
        lambda_i <- exp(drop(X[i, , drop = FALSE] %*% BetaMat))
        S[, , i] <- Q %*% diag(lambda_i, p) %*% t(Q)
    }

    list(
        S = S,
        X = X,
        Q = Q,
        BetaMat = BetaMat,
        H = apply(S, c(1, 2), mean),
        p = p,
        n = n
    )
}
