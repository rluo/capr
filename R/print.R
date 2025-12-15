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
#'
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
        cat("Log-likelihood by component: ",
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

    cat("\nCoefficient estimates (B-hat):\n")
    coef_df <- data.frame(
        Term = row_labs,
        coef_fmt,
        check.names = FALSE,
        row.names = NULL
    )
    print(coef_df, row.names = FALSE, right = TRUE, ...)
    invisible(x)
}
