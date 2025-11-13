test_that("capr estimates multiple orthogonal components", {
    sim <- capr_test_data(seed = 123L, p = 5L, n = 120L)
    K <- 2L

    fit <- capr(
        S = sim$S,
        X = sim$X,
        K = K,
        max_iter = 250L,
        tol = 1e-8,
        orth = TRUE
    )

    expect_equal(dim(fit$B), c(ncol(sim$X), K))
    expect_equal(dim(fit$Gamma), c(sim$p, K))
    expect_true(all(is.finite(fit$B)))
    expect_true(all(is.finite(fit$Gamma)))

    H <- sim$H
    gram <- t(fit$Gamma) %*% H %*% fit$Gamma
    expect_equal(diag(gram), rep(1, K), tolerance = 1e-6)
    expect_lt(abs(gram[1, 2]), 1e-6)

    align_scores <- cosine_alignment_scores(fit$Gamma, sim$Q)
    expect_true(all(align_scores > 0.75))
})
