test_that("capr::single component::no init::orth", {
    simu.data <- simu.capr(seed = 123L, p = 5L, n = 120L)


    fit <- capr(
        S = simu.data$S,
        X = simu.data$X,
        K = 1L,
        max_iter = 200L,
        tol = 1e-6,
        orth = TRUE
    )

    expect_equal(dim(fit$B), c(ncol(simu.data$X), 1L))
    expect_equal(dim(fit$Gamma), c(simu.data$p, 1L))
    expect_true(all(is.finite(fit$B)))
    expect_true(all(is.finite(fit$Gamma)))

    ## normalization check
    gamma_hat <- fit$Gamma[, 1]
    norm_check <- drop(t(gamma_hat) %*% simu.data$H %*% gamma_hat)
    expect_equal(norm_check, 1, tolerance = 1e-6)

    ## estimate gamma
    align_score <- cosine_alignment_scores(fit$Gamma, simu.data$Q)
    expect_gt(align_score[1], 0.8)

    ## estimate beta
    expect_lt(abs(fit$B[2, 1] - simu.data$B[2, 2]), 0.3) # true beta value
})



test_that("capr::single component::no init::no orth", {
    simu.data <- simu.capr(seed = 123L, p = 5L, n = 120L)

    fit <- capr(
        S = simu.data$S,
        X = simu.data$X,
        K = 1L,
        max_iter = 200L,
        tol = 1e-6,
        orth = FALSE
    )

    expect_equal(dim(fit$B), c(ncol(simu.data$X), 1L))
    expect_equal(dim(fit$Gamma), c(simu.data$p, 1L))
    expect_true(all(is.finite(fit$B)))
    expect_true(all(is.finite(fit$Gamma)))

    ## normalization check
    gamma_hat <- fit$Gamma[, 1]
    norm_check <- drop(t(gamma_hat) %*% simu.data$H %*% gamma_hat)
    expect_equal(norm_check, 1, tolerance = 1e-6)

    ## estimate gamma
    align_score <- cosine_alignment_scores(fit$Gamma, simu.data$Q)
    expect_gt(align_score[1], 0.8)

    ## estimate beta
    expect_lt(abs(fit$B[2, 1] - simu.data$B[2, 2]), 0.3) # true beta value
})



test_that("capr::single component::init::orth", {
    simu.data <- simu.capr(seed = 123L, p = 5L, n = 120L)


    fit <- capr(
        S = simu.data$S,
        X = simu.data$X,
        K = 1L,
        B.init = array(c(0, simu.data$B[2, 2]), dim = c(ncol(simu.data$X), 10L, 1L)),
        Gamma.init = array(simu.data$Q[, 1, drop = TRUE], dim = c(simu.data$p, 10L, 1L)),
        max_iter = 200L,
        tol = 1e-6,
        orth = TRUE,
        n.init = 10L
    )

    expect_equal(dim(fit$B), c(ncol(simu.data$X), 1L))
    expect_equal(dim(fit$Gamma), c(simu.data$p, 1L))
    expect_true(all(is.finite(fit$B)))
    expect_true(all(is.finite(fit$Gamma)))

    ## normalization check
    gamma_hat <- fit$Gamma[, 1]
    norm_check <- drop(t(gamma_hat) %*% simu.data$H %*% gamma_hat)
    expect_equal(norm_check, 1, tolerance = 1e-6)

    ## estimate gamma
    align_score <- cosine_alignment_scores(fit$Gamma, simu.data$Q)
    expect_gt(align_score[1], 0.8)

    ## estimate beta
    expect_lt(abs(fit$B[2, 1] - simu.data$B[2, 2]), 0.3) # true beta value
})



test_that("capr::single component::init::no orth", {
    simu.data <- simu.capr(seed = 123L, p = 5L, n = 120L)


    fit <- capr(
        S = simu.data$S,
        X = simu.data$X,
        K = 1L,
        B.init = array(c(0, simu.data$B[2, 2]), dim = c(ncol(simu.data$X), 10L, 1L)),
        Gamma.init = array(simu.data$Q[, 1, drop = TRUE], dim = c(simu.data$p, 10L, 1L)),
        max_iter = 200L,
        tol = 1e-6,
        orth = FALSE,
        n.init = 10L
    )

    expect_equal(dim(fit$B), c(ncol(simu.data$X), 1L))
    expect_equal(dim(fit$Gamma), c(simu.data$p, 1L))
    expect_true(all(is.finite(fit$B)))
    expect_true(all(is.finite(fit$Gamma)))

    ## normalization check
    gamma_hat <- fit$Gamma[, 1]
    norm_check <- drop(t(gamma_hat) %*% simu.data$H %*% gamma_hat)
    expect_equal(norm_check, 1, tolerance = 1e-6)

    ## estimate gamma
    align_score <- cosine_alignment_scores(fit$Gamma, simu.data$Q)
    expect_gt(align_score[1], 0.8)

    ## estimate beta
    expect_lt(abs(fit$B[2, 1] - simu.data$B[2, 2]), 0.3) # true beta value
})
