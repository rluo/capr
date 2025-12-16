test_that("capr_bootstrap produces coefficient intervals", {
    sim <- capr_test_data(seed = 321L, p = 5L, n = 80L)
    K <- 2L

    fit <- capr(
        S = sim$S,
        X = sim$X,
        K = K,
        max_iter = 150L,
        tol = 1e-7,
        orth = TRUE
    )

    boot <- capr.boot(
        S = sim$S,
        X = sim$X,
        fit = fit,
        nboot = 6L,
        max_iter = 20L,
        tol = 1e-6,
        seed = 123L
    )

    expect_true(inherits(boot, "capr.boot"))
    expect_equal(dim(boot$beta), c(ncol(sim$X), K))
    expect_equal(dim(boot$ci_lower), c(ncol(sim$X), K))
    expect_equal(dim(boot$ci_upper), c(ncol(sim$X), K))
    expect_equal(boot$level, 0.95)
    expect_true(all(is.finite(boot$beta)))
    expect_true(all(boot$ci_lower <= boot$ci_upper))
})
