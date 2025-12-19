#' Simulate covariance matrices compatible with `capr()`
#'
#' Generates a simple synthetic dataset for CAP regression consisting of a
#' covariance cube, design matrix, and the latent orthogonal directions used to
#' build the covariance slices.
#'
#' @param seed Integer seed used for reproducibility.
#' @param p Number of variables (dimension of the covariance matrices).
#' @param n Number of observations (slices) to generate.
#'
#' @return A list with components:
#' \item{S}{Array of dimension \eqn{p \times p \times n} holding the simulated
#'   covariance matrices.}
#' \item{X}{Design matrix of size \eqn{n \times 2} with an intercept and a
#'   Bernoulli covariate.}
#' \item{Q}{Orthogonal matrix whose columns are the latent directions.}
#' \item{BetaMat}{True coefficient matrix used to form the eigenvalues.}
#' \item{H}{Average covariance matrix \eqn{\frac{1}{n}\sum_i S_i}.}
#' \item{p, n}{The dimension and sample size supplied to the generator.}
#'
#' @examples
#' sim <- simu.capr(seed = 10, p = 4, n = 50)
#' str(sim$S)
#'
#' @export
simu.capr <- function(seed = 123L, p = 5L, n = 120L) {
    stopifnot(length(seed) == 1L, length(p) == 1L, length(n) == 1L)
    set.seed(as.integer(seed))

    BetaMat <- rbind(
        c(5, 4, 1, -1, -2),
        c(0, -1.2, 0.9, 0, 0)
    )
    if (ncol(BetaMat) != p) {
        stop("BetaMat must have ", p, " columns.", call. = FALSE)
    }

    Q <- qr.Q(qr(matrix(stats::rnorm(p * p), p, p)))
    X <- cbind(1, stats::rbinom(n, size = 1L, prob = 0.5))
    colnames(X) <- c("Intercept", "Covariate")
    S <- array(0, dim = c(p, p, n))
    for (i in seq_len(n)) {
        lambda_i <- exp(drop(X[i, , drop = FALSE] %*% BetaMat))
        Sigma <- Q %*% diag(lambda_i, p) %*% t(Q)
        Y <- MASS::mvrnorm(n = 200L, mu = rep(0, p), Sigma = Sigma)
        S[, , i] <- stats::cov(Y)
        ## S[, , i] <- Q %*% diag(lambda_i, p) %*% t(Q)
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
