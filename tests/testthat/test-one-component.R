test_that("capr returns a normalized direction for a single component", {
    sim <- capr_test_data(seed = 42L, p = 5L, n = 100L)

    fit <- capr(
        S = sim$S,
        X = sim$X,
        K = 1L,
        max_iter = 200L,
        tol = 1e-8,
        orth = TRUE
    )

    expect_equal(dim(fit$B), c(ncol(sim$X), 1L))
    expect_equal(dim(fit$Gamma), c(sim$p, 1L))
    expect_true(all(is.finite(fit$B)))
    expect_true(all(is.finite(fit$Gamma)))

    gamma_hat <- fit$Gamma[, 1]
    norm_check <- drop(t(gamma_hat) %*% sim$H %*% gamma_hat)
    expect_equal(norm_check, 1, tolerance = 1e-6)

    align_score <- cosine_alignment_scores(fit$Gamma, sim$Q)
    expect_gt(align_score[1], 0.8)
})
