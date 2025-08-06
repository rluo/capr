#' First CAP projection direction (Algorithm 1)
#'
#' @param Y  list of numeric matrices; each matrix is T_i × p for subject i
#' @param X  n × q covariate matrix (first column must be 1’s for intercept)
#' @param tol        convergence tolerance (default 1e-6)
#' @param max_iter   maximum iterations (default 1000)
#' @param beta_init  optional length-q numeric vector
#' @param gamma_init optional length-p numeric vector
#'
#' @return list with elements beta, gamma, objective, iterations, converged
#' @examples
#' ## toy example
#' set.seed(1)
#' n <- 10
#' p <- 4
#' Ti <- 30
#' X <- cbind(1, rnorm(n))
#' Y <- lapply(seq_len(n), function(i) matrix(rnorm(Ti * p), Ti, p))
#' res <- cap_first_direction(Y, X)
#' res$beta
#' res$gamma
#' @export
cap_first_direction <- function(Y, X,
                                tol = 1e-6, max_iter = 1000,
                                beta_init = NULL,
                                gamma_init = NULL) {
    stopifnot(is.list(Y), is.matrix(X))
    n <- length(Y)
    p <- ncol(Y[[1]])
    q <- ncol(X)
    if (is.null(beta_init)) beta_init <- rep(0, q)
    if (is.null(gamma_init)) {
        # leading eigenvector of pooled covariance – reasonable default
        Ssum <- Reduce(`+`, lapply(Y, crossprod))
        gamma_init <- eigen(Ssum)$vectors[, 1]
    }
    cap_first_direction_cpp(
        Y_list = Y,
        X = X,
        beta = beta_init,
        gamma = gamma_init,
        tol = tol,
        max_iter = max_iter
    )
}
