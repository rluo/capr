project_component_variance <- function(S_cube, gamma) {
    vapply(
        seq_len(dim(S_cube)[3L]),
        function(i) {
            as.numeric(crossprod(gamma, S_cube[, , i] %*% gamma))
        },
        numeric(1L)
    )
}

solve_beta_fixed_gamma <- function(y, X, weight, beta_start, max_iter, tol) {
    beta <- beta_start
    for (iter in seq_len(max_iter)) {
        eta <- as.vector(X %*% beta)
        wi <- exp(-eta) * y
        wi[!is.finite(wi) | wi <= 0] <- .Machine$double.eps
        g <- crossprod(X, weight - wi)
        WX <- sqrt(wi) * X
        H <- crossprod(WX)
        step <- solve_symmetric(H, g)
        beta_new <- beta - step
        if (max(abs(beta_new - beta)) < tol) {
            beta <- beta_new
            break
        }
        beta <- beta_new
    }
    beta
}

solve_symmetric <- function(H, g) {
    chol_attempt <- tryCatch(chol(H), error = function(e) NULL)
    if (!is.null(chol_attempt)) {
        return(backsolve(chol_attempt, forwardsolve(t(chol_attempt), g)))
    }
    qr.solve(H, g)
}


#' Bootstrap confidence intervals for CAP coefficients
#'
#' Generates bootstrap replicates of the CAP regression coefficients while
#' holding the fitted directions \eqn{\Gamma} fixed. Each replicate samples the
#' covariance slices \eqn{S_i} with replacement, projects them onto the fixed
#' directions to obtain component-specific variances, and re-solves the
#' \eqn{\beta^{(k)}} equations. Quantile-based confidence intervals are
#' returned for every predictor/component pair.
#'
#' @param S Numeric array of shape \eqn{p \times p \times n}, matching the data
#'   used to fit \code{fit}.
#' @param X Numeric design matrix with \eqn{n} rows.
#' @param fit A \code{\link{capr}} fit containing \code{$B} and \code{$Gamma}.
#' @param nboot Number of bootstrap replicates.
#' @param weight Optional numeric vector of length \eqn{n} (defaults to the
#'   weights used during the original fit, or 1).
#' @param level Confidence level for the returned intervals.
#' @param max_iter Maximum Newton iterations for solving \eqn{\beta}.
#' @param tol Convergence tolerance for the Newton solver.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with:
#' \item{samples}{Array of dimension \eqn{q \times K \times nboot} storing the
#'   bootstrap coefficient estimates.}
#' \item{ci_lower, ci_upper}{Matrices \eqn{q \times K} with the lower/upper
#'   confidence limits.}
#' \item{level}{The requested confidence level.}
#'
#' @export
capr_bootstrap <- function(S, X, fit, nboot = 200L, weight = NULL,
                           level = 0.95, max_iter = 50L, tol = 1e-6,
                           seed = NULL) {
    if (!inherits(fit, "capr")) {
        stop("`fit` must be an object returned by `capr()`.", call. = FALSE)
    }
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

    if (is.null(weight)) {
        weight <- rep(1, n)
    } else if (length(weight) != n) {
        stop("`weight` must have length n.", call. = FALSE)
    }
    if (!is.numeric(weight)) {
        stop("`weight` must be numeric.", call. = FALSE)
    }

    samples <- array(NA_real_, dim = c(q, K, nboot),
                     dimnames = list(rownames(fit$B), colnames(fit$B), NULL))

    for (b in seq_len(nboot)) {
        idx <- sample.int(n, n, replace = TRUE)
        S_boot <- S[, , idx, drop = FALSE]
        X_boot <- X[idx, , drop = FALSE]
        w_boot <- weight[idx]

        for (k in seq_len(K)) {
            gamma_k <- Gamma_hat[, k]
            beta_start <- fit$B[, k]
            y_proj <- project_component_variance(S_boot, gamma_k)
            beta_hat <- solve_beta_fixed_gamma(
                y = y_proj,
                X = X_boot,
                weight = w_boot,
                beta_start = beta_start,
                max_iter = max_iter,
                tol = tol
            )
            samples[, k, b] <- beta_hat
        }
    }

    alpha <- 1 - level
    lower <- apply(samples, c(1, 2), stats::quantile, probs = alpha / 2, na.rm = TRUE)
    upper <- apply(samples, c(1, 2), stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)

    list(samples = samples, ci_lower = lower, ci_upper = upper, level = level)
}
