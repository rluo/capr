test_that("capr_bootstrap produces coefficient intervals", {
    sim <- capr_test_data(seed = 11L, n = 60L)
    fit <- capr(
        S = sim$S,
        X = sim$X,
        K = 1L,
        n.init = 3L,
        max_iter = 50L,
        tol = 1e-6
    )

    boot <- capr_bootstrap(
        S = sim$S,
        X = sim$X,
        fit = fit,
        nboot = 6L,
        max_iter = 10L,
        tol = 1e-6,
        seed = 99L
    )

    expect_equal(dim(boot$samples), c(ncol(sim$X), 1L, 6L))
    expect_equal(dim(boot$ci_lower), c(ncol(sim$X), 1L))
    expect_equal(dim(boot$ci_upper), c(ncol(sim$X), 1L))
    expect_true(all(boot$ci_lower <= boot$ci_upper))
})
