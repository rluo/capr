#' Plot deviation diagnostics by component count
#'
#' For a fitted CAP regression, plots two diagnostics across the first
#' \eqn{K} components: (1) the negative log-likelihood returned by [capr()]
#' and (2) the log deviation-from-diagonality (DfD) for the loading matrix
#' formed by the first \eqn{k} directions. Both curves help assess the gain
#' from adding components.
#'
#' @param x A `capr` object returned by [capr()].
#' @param ... Additional arguments passed to [graphics::plot()] and applied to
#'   both panels (for example, `pch`, `col`, or axis limits).
#'
#' @details
#' The DfD criterion for the first \eqn{k} directions \eqn{\Gamma^{(k)}} is
#' \deqn{\text{DfD}(\Gamma^{(k)}) =
#' \left(\prod_{i = 1}^{n}
#'   \nu\!\left(\Gamma^{(k)\top} S_i \Gamma^{(k)} / T_i\right)^{T_i}
#' \right)^{1 / \sum_i T_i},}
#' where
#' \deqn{
#' \nu(A)=\frac{\det\{\mathrm{diag}(A)\}}{\det(A)}
#' } for a positive definite matrix \eqn{A}.
#' The curve shows \eqn{\log \text{DfD}(\Gamma^{(k)})}. A common choice for
#' \eqn{k} is the last point before a sudden jump in the negative
#' log-likelihood or log-DfD curve.
#'
#' @return Invisibly returns the numeric vector of log deviation values (one per
#'   component).
#'
#' @seealso [log_deviation_from_diagonality()]
#' @examples
#' \dontrun{
#' sim <- simu.capr(seed = 123L, n = 120L)
#' fit <- capr(S = sim$S, X = sim$X, K = 3L)
#' plot(fit)
#' }
#' @method plot capr
#' @export
plot.capr <- function(x, ...) {
    Gamma <- x$Gamma
    K <- ncol(x$Gamma)
    weight <- x$weight
    logdfd <- rep(0, K)
    S <- x$S
    for (k in seq_len(K)) {
        Gamma_k <- Gamma[, 1:k, drop = FALSE]
        logdfd[k] <- log_deviation_from_diagonality(S, weight, Gamma_k)
    }
    # Save current graphics parameters to restore later
    op <- graphics::par(ask = TRUE)
    # Ensure parameters are reset when the function exits (good practice)
    on.exit(graphics::par(op))
    plot(
        x = seq_along(logdfd),
        y = x$loglike,
        type = "b",
        xlab = "Number of Components",
        ylab = "CAP Neg Log-Likelihood",
        main = "CAP Neg Log-Likelihood by Number of Components",
        ...
    )
    plot(
        x = seq_along(logdfd),
        y = logdfd,
        type = "b",
        xlab = "Number of Components",
        ylab = "Log-Deviance from Diagonality",
        main = "CAP Log-Deviance from Diagonality",
        ...
    )

    graphics::par(op) # Restore original graphics parameters
    invisible(logdfd)
}
