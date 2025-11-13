# -----------------------------------------------------------------------------
#  Flury-Gautschi (1986) example reproduced via FG() + log deviation criterion
# -----------------------------------------------------------------------------

test_that("FG decreases the log-deviation criterion in the classic example", {
    skip_on_cran()

    P <- 6L
    M <- 2L
    cov_cube <- array(0, dim = c(P, P, M))

    cov_cube[, , 1] <- matrix(c(
        45, 10, 0, 5, 0, 0,
        10, 45, 5, 0, 0, 0,
        0, 5, 45, 10, 0, 0,
        5, 0, 10, 45, 0, 0,
        0, 0, 0, 0, 16.4, -4.8,
        0, 0, 0, 0, -4.8, 13.6
    ), P, P, byrow = TRUE)

    cov_cube[, , 2] <- matrix(c(
        27.5, -12.5, -0.5, -4.5, -2.04, 3.72,
        -12.5, 27.5, -4.5, -0.5, 2.04, -3.72,
        -0.5, -4.5, 24.5, -9.5, -3.72, -2.04,
        -4.5, -0.5, -9.5, 24.5, 3.72, 2.04,
        -2.04, 2.04, -3.72, 3.72, 54.76, -4.68,
        3.72, -3.72, -2.04, 2.04, -4.68, 51.24
    ), P, P, byrow = TRUE)

    B_ref <- matrix(c(
        0.5, -0.5545, 0.5, 0.4281, -0.0956, 0.0083,
        0.5, 0.5545, 0.5, -0.4281, 0.0956, -0.0083,
        -0.5, -0.4250, 0.5, -0.5623, -0.0541, -0.0169,
        -0.5, 0.4250, 0.5, 0.5623, 0.0541, 0.0169,
        0.0, -0.1265, 0.0, 0.0014, 0.7919, 0.5974,
        0.0, 0.0878, 0.0, -0.0337, -0.5906, 0.8015
    ), P, P, byrow = TRUE)

    fit <- FG(cov_cube, maxit = 200L)

    expect_equal(dim(fit), c(P, P))
    expect_true(all(is.finite(fit)))
    expect_equal(crossprod(fit), diag(P), tolerance = 1e-6)

    weights <- rep(1, M)
    crit_identity <- log_deviation_from_diagonality(cov_cube, weights, diag(P))
    crit_fit <- log_deviation_from_diagonality(cov_cube, weights, fit)
    crit_ref <- log_deviation_from_diagonality(cov_cube, weights, B_ref)

    expect_equal(crit_ref, 0.03432332, tolerance = 1e-6)
    expect_lt(crit_fit, crit_identity)
    expect_lt(crit_fit, crit_identity * 0.5)
    expect_lt(crit_fit, crit_ref + 0.25)
})
