#' Bootstrap confidence intervals for CAP coefficients
#'
#' Generates bootstrap inference  of the CAP regression coefficients while
#' holding the fitted directions \eqn{\Gamma} fixed. Each replicate samples the
#' covariance slices \eqn{S[,,i]} with replacement, projects them onto the fixed
#' directions to obtain component-specific variances, and re-solves the
#' \eqn{\beta^{(k)}} equations. Quantile-based confidence intervals are
#' returned for every predictor/component pair.
#'
#' @param fit A \code{\link{capr}} fit containing \code{$B} and \code{$Gamma}.
#' @param nboot Number of bootstrap replicates.
#' @param level Confidence level for the returned intervals.
#' @param max_iter Maximum Newton iterations for solving \eqn{\beta}.
#' @param tol Convergence tolerance for the Newton solver.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with:
#' \item{beta}{bootstrap average of beta of dimension \eqn{q \times K} }
#' \item{ci_lower, ci_upper}{Matrices \eqn{q \times K} with the lower/upper
#'   confidence limits.}
#' \item{level}{The requested confidence level.}
#'
#' @export
capr.boot <- function(fit, nboot = 1000L,
                      level = 0.95, max_iter = 100L, tol = 1e-6,
                      seed = NULL) {
    if (!inherits(fit, "capr")) {
        stop("`fit` must be an object returned by `capr()`.", call. = FALSE)
    }

    S <- fit$S
    X <- fit$X

    if (!is.array(S) || length(dim(S)) != 3L) {
        stop("`S` must be a numeric 3D array.", call. = FALSE)
    }
    if (!is.matrix(X) || !is.numeric(X)) {
        stop("`X` must be a numeric matrix.", call. = FALSE)
    }
    if (nrow(X) != dim(S)[3L]) {
        stop("nrow(X) must equal dim(S)[3].", call. = FALSE)
    }

    n <- dim(S)[3L]
    p <- dim(S)[1L]
    q <- ncol(X)
    K <- ncol(fit$B)

    if (!is.null(seed)) set.seed(as.integer(seed))

    Gamma_hat <- fit$Gamma
    if (!is.matrix(Gamma_hat) || nrow(Gamma_hat) != p || ncol(Gamma_hat) != K) {
        stop("`fit$Gamma` must be a p x K matrix compatible with `S`.", call. = FALSE)
    }
    if (!is.matrix(fit$B) || nrow(fit$B) != q || ncol(fit$B) != K) {
        stop("`fit$B` must be a q x K matrix compatible with `X`.", call. = FALSE)
    }

    weight <- fit$weight


    beta_all <- array(0, dim = c(q, K, nboot))

    for (b in seq_len(nboot)) {
        idx <- sample.int(n, n, replace = TRUE)

        for (k in seq_len(K)) {
            X_boot <- X[idx, , drop = FALSE]
            w_boot <- weight[idx]
            beta_start <- fit$B[, k]
            beta_hat <- newton_beta(
                S = S[, , idx],
                X = X_boot,
                T = w_boot,
                beta_init = beta_start,
                gamma = Gamma_hat[, k],
                max_iter = max_iter,
                tol = tol
            )
            beta_all[, k, b] <- beta_hat
        }
    }

    alpha <- 1 - level
    lower <- apply(beta_all, c(1, 2), stats::quantile, probs = alpha / 2, na.rm = TRUE)
    upper <- apply(beta_all, c(1, 2), stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)

    ret <- list(beta = apply(beta_all, c(1, 2), mean, na.rm = TRUE), ci_lower = lower, ci_upper = upper, level = level)
    class(ret) <- c("capr.boot")
    return(ret)
}
