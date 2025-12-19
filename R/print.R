#' Print method for CAP regression fits
#'
#' Formats the coefficient matrix \eqn{\hat{B}} returned by [capr()] in a
#' linear-regression style table, showing the estimate for each predictor and
#' component.
#'
#' @param x An object of class `capr`, typically the result of [capr()].
#' @param digits Number of significant digits to show when printing numeric
#'   values.
#' @param ... Additional arguments passed on to [print.data.frame()].
#'
#' @return The input object `x`, invisibly.
#' @examples
#' simu.data <- simu.capr(seed = 123L, p = 5L, n = 120L)
#' K <- 2L
#' fit <- capr(
#'     S = simu.data$S,
#'     X = simu.data$X,
#'     K = K
#' )
#' print(fit)
#' @export
print.capr <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    if (!inherits(x, "capr")) {
        stop("`x` must be a capr fit produced by `capr()`.", call. = FALSE)
    }
    B_hat <- x$B
    Gamma_hat <- x$Gamma
    if (is.null(B_hat) || is.null(Gamma_hat)) {
        cat("capr object with no coefficient or direction estimates.\n")
        return(invisible(x))
    }

    if (!is.matrix(B_hat)) {
        B_hat <- as.matrix(B_hat)
    }
    if (!is.matrix(Gamma_hat)) {
        Gamma_hat <- as.matrix(Gamma_hat)
    }

    K <- ncol(B_hat)
    cat(sprintf(
        "Covariate Assisted Projection Regression with %d component%s.\n",
        K,
        if (K == 1L) "" else "s"
    ))

    if (!is.null(x$loglike)) {
        loglike_fmt <- formatC(
            x$loglike,
            digits = digits,
            format = "f"
        )
        cat("Objective function by component: ",
            paste(loglike_fmt, collapse = ", "),
            "\n",
            sep = ""
        )
    }

    row_labs <- rownames(B_hat)
    if (is.null(row_labs)) {
        row_labs <- paste0("x", seq_len(nrow(B_hat)))
    }
    col_labs <- colnames(B_hat)
    if (is.null(col_labs)) {
        col_labs <- paste0("Comp", seq_len(ncol(B_hat)))
    }

    coef_fmt <- signif(B_hat, digits = digits)
    dimnames(coef_fmt) <- list(row_labs, col_labs)

    cat("\nCoefficient estimates:\n")
    coef_df <- data.frame(
        Term = row_labs,
        coef_fmt,
        check.names = FALSE,
        row.names = NULL
    )
    print(coef_df, row.names = FALSE, right = TRUE, ...)
    invisible(x)
}

#' Print method for `capr.boot` objects
#'
#' Displays bootstrap coefficient estimates and their confidence intervals
#' component by component as compact tables.
#'
#' @param x An object of class `capr.boot`, typically produced by [capr.boot()].
#' @inheritParams print.capr
#'
#' @return The input object `x`, invisibly.
#' @examples
#' simu.data <- simu.capr(seed = 123L, p = 5L, n = 120L)
#' K <- 2L
#' fit <- capr(
#'     S = simu.data$S,
#'     X = simu.data$X,
#'     K = K
#' )
#' fit.boot <- capr.boot(
#'     S = simu.data$S,
#'     X = simu.data$X,
#'     fit = fit,
#'     nboot = 10L,
#'     max_iter = 20L,
#'     tol = 1e-6,
#'     seed = 123L
#' )
#' print(fit.boot)
#' @export
print.capr.boot <- function(x, digits = max(4L, getOption("digits") - 4L), ...) {
    if (!inherits(x, "capr.boot")) {
        stop("`x` must be a bootstrap result from `capr.boot()`.", call. = FALSE)
    }
    beta_hat <- x$beta
    ci_lower <- x$ci_lower
    ci_upper <- x$ci_upper
    if (is.null(beta_hat) || is.null(ci_lower) || is.null(ci_upper)) {
        cat("capr.boot object with no coefficient summaries.\n")
        return(invisible(x))
    }

    if (!is.matrix(beta_hat)) beta_hat <- as.matrix(beta_hat)
    if (!is.matrix(ci_lower)) ci_lower <- as.matrix(ci_lower)
    if (!is.matrix(ci_upper)) ci_upper <- as.matrix(ci_upper)

    K <- ncol(beta_hat)
    q <- nrow(beta_hat)
    if (!nrow(ci_lower) == q || !ncol(ci_lower) == K ||
        !nrow(ci_upper) == q || !ncol(ci_upper) == K) {
        stop("Bootstrap coefficient and interval matrices must have matching dimensions.", call. = FALSE)
    }

    row_labs <- rownames(beta_hat)
    if (is.null(row_labs)) {
        row_labs <- paste0("x", seq_len(q))
    }
    col_labs <- colnames(beta_hat)
    if (is.null(col_labs)) {
        col_labs <- paste0("Comp", seq_len(K))
    }

    level <- x$level
    if (is.null(level) || !is.finite(level)) {
        level <- NA_real_
    }

    level_text <- if (is.na(level)) {
        "Bootstrap coefficient summaries"
    } else {
        sprintf("CAPR coefficient summaries (level %.1f%%)", 100 * level)
    }
    cat(level_text, "\n")

    for (k in seq_len(K)) {
        cat("\n", col_labs[k], ":\n", sep = "")
        comp_df <- data.frame(
            Term = row_labs,
            Estimate = signif(beta_hat[, k], digits = digits),
            CI_lower = signif(ci_lower[, k], digits = digits),
            CI_upper = signif(ci_upper[, k], digits = digits),
            row.names = NULL,
            check.names = FALSE
        )
        print(comp_df, row.names = FALSE, right = TRUE, ...)
    }
    invisible(x)
}
