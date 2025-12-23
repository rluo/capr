test_that("capr::multiple components:no init::orth", {
    sim <- simu.capr(seed = 123L, p = 5L, n = 150L)

    K <- 2L

    fit <- capr(
        S = sim$S,
        X = sim$X,
        K = K,
        max_iter = 250L,
        tol = 1e-4,
        orth = TRUE
    )

    expect_equal(dim(fit$B), c(ncol(sim$X), K))
    expect_equal(dim(fit$Gamma), c(sim$p, K))
    expect_true(all(is.finite(fit$B)))
    expect_true(all(is.finite(fit$Gamma)))

    ## orthogonality check
    H <- sim$H
    gram <- t(fit$Gamma) %*% H %*% fit$Gamma
    expect_equal(unname(diag(gram)), rep(1, K), tolerance = 1e-3)
    expect_lt(abs(gram[1, 2]), 1e-3)

    ## estimate gamma
    align_scores <- cosine_alignment_scores(fit$Gamma, sim$Q)
    expect_true(all(align_scores > 0.75))

    ## estimate beta
    expect_lt(abs(fit$B[2, 1] - sim$BetaMat[2, 2]), 0.3) # true beta value
    expect_lt(abs(fit$B[2, 2] - sim$BetaMat[2, 3]), 0.3) # true beta value
})



test_that("capr::multiple components:init::orth", {
    sim <- simu.capr(seed = 123L, p = 5L, n = 150L)

    K <- 2L

    B.init <- array(0, dim = c(ncol(sim$X), 10L, K))
    B.init[, , 1] <- matrix(c(0, sim$BetaMat[2, 2]), nrow = ncol(sim$X), ncol = 10L)
    B.init[, , 2] <- matrix(c(0, sim$BetaMat[2, 3]), nrow = ncol(sim$X), ncol = 10L)
    Gamma.init <- array(0, dim = c(sim$p, 10L, K))
    Gamma.init[, , 1] <- matrix(sim$Q[, 2, drop = TRUE], nrow = sim$p, ncol = 10L)
    Gamma.init[, , 2] <- matrix(sim$Q[, 3, drop = TRUE], nrow = sim$p, ncol = 10L)
    fit <- capr(
        S = sim$S,
        X = sim$X,
        K = K,
        B.init = B.init,
        Gamma.init = Gamma.init,
        max_iter = 250L,
        tol = 1e-4,
        orth = TRUE,
        n.init = 10L
    )

    expect_equal(dim(fit$B), c(ncol(sim$X), K))
    expect_equal(dim(fit$Gamma), c(sim$p, K))
    expect_true(all(is.finite(fit$B)))
    expect_true(all(is.finite(fit$Gamma)))

    ## orthogonality check
    H <- sim$H
    gram <- t(fit$Gamma) %*% H %*% fit$Gamma
    expect_equal(unname(diag(gram)), rep(1, K), tolerance = 1e-3)
    expect_lt(abs(gram[1, 2]), 1e-3)

    ## estimate gamma
    align_scores <- cosine_alignment_scores(fit$Gamma, sim$Q)
    expect_true(all(align_scores > 0.75))

    ## estimate beta
    expect_lt(abs(fit$B[2, 1] - sim$BetaMat[2, 2]), 0.3) # true beta value
    expect_lt(abs(fit$B[2, 2] - sim$BetaMat[2, 3]), 0.3) # true beta value
})




# test_that("capr::multiple components:no init::no orth", {
#     sim <- simu.capr(seed = 123L, p = 5L, n = 150L)

#     K <- 2L

#     fit <- capr(
#         S = sim$S,
#         X = sim$X,
#         K = K,
#         max_iter = 250L,
#         tol = 1e-4,
#         orth = FALSE
#     )

#     expect_equal(dim(fit$B), c(ncol(sim$X), K))
#     expect_equal(dim(fit$Gamma), c(sim$p, K))
#     expect_true(all(is.finite(fit$B)))
#     expect_true(all(is.finite(fit$Gamma)))

#     ## norm check
#     H <- sim$H
#     gram <- t(fit$Gamma) %*% H %*% fit$Gamma
#     expect_equal(unname(diag(gram)), rep(1, K), tolerance = 1e-1)
#     # expect_lt(abs(gram[1, 2]), 1e-3)

#     ## estimate gamma
#     align_scores <- cosine_alignment_scores(fit$Gamma, sim$Q)
#     expect_true(all(align_scores > 0.75))

#     ## estimate beta
#     expect_lt(abs(fit$B[2, 1] - sim$BetaMat[2, 2]), 0.3) # true beta value
#     expect_lt(abs(fit$B[2, 2] - sim$BetaMat[2, 3]), 0.3) # true beta value
# })



test_that("capr::multiple components:init::no orth", {
    sim <- simu.capr(seed = 123L, p = 5L, n = 150L)

    K <- 2L

    B.init <- array(0, dim = c(ncol(sim$X), 10L, K))
    B.init[, , 1] <- matrix(c(0, sim$BetaMat[2, 2]), nrow = ncol(sim$X), ncol = 10L)
    B.init[, , 2] <- matrix(c(0, sim$BetaMat[2, 3]), nrow = ncol(sim$X), ncol = 10L)
    Gamma.init <- array(0, dim = c(sim$p, 10L, K))
    Gamma.init[, , 1] <- matrix(sim$Q[, 2, drop = TRUE], nrow = sim$p, ncol = 10L)
    Gamma.init[, , 2] <- matrix(sim$Q[, 3, drop = TRUE], nrow = sim$p, ncol = 10L)
    fit <- capr(
        S = sim$S,
        X = sim$X,
        K = K,
        B.init = B.init,
        Gamma.init = Gamma.init,
        max_iter = 250L,
        tol = 1e-4,
        orth = FALSE,
        n.init = 10L
    )

    expect_equal(dim(fit$B), c(ncol(sim$X), K))
    expect_equal(dim(fit$Gamma), c(sim$p, K))
    expect_true(all(is.finite(fit$B)))
    expect_true(all(is.finite(fit$Gamma)))

    ## norm check
    H <- sim$H
    gram <- t(fit$Gamma) %*% H %*% fit$Gamma
    expect_equal(unname(diag(gram)), rep(1, K), tolerance = 1e-3)
    # expect_lt(abs(gram[1, 2]), 1e-3)

    ## estimate gamma
    align_scores <- cosine_alignment_scores(fit$Gamma, sim$Q)
    expect_true(all(align_scores > 0.75))

    ## estimate beta
    expect_lt(abs(fit$B[2, 1] - sim$BetaMat[2, 2]), 0.3) # true beta value
    expect_lt(abs(fit$B[2, 2] - sim$BetaMat[2, 3]), 0.3) # true beta value
})
