capr_test_data <- function(seed = 123L, p = 5L, n = 120L) {
    stopifnot(length(seed) == 1L, length(p) == 1L, length(n) == 1L)
    set.seed(as.integer(seed))

    BetaMat <- rbind(
        c(5, 4, 1, -1, -2),
        c(0, -1, 1, 0, 0)
    )
    if (ncol(BetaMat) != p) {
        stop("BetaMat must have ", p, " columns.", call. = FALSE)
    }

    Q <- qr.Q(qr(matrix(stats::rnorm(p * p), p, p)))
    X <- cbind(1, stats::rbinom(n, size = 1L, prob = 0.5))

    S <- array(0, dim = c(p, p, n))
    for (i in seq_len(n)) {
        lambda_i <- exp(drop(X[i, , drop = FALSE] %*% BetaMat))
        S[, , i] <- Q %*% diag(lambda_i, p) %*% t(Q)
    }

    list(
        S = S,
        X = X,
        Q = Q,
        BetaMat = BetaMat,
        H = apply(S, c(1, 2), mean),
        p = p,
        n = n
    )
}

cosine_similarity_r <- function(a, b) {
    num <- sum(a * b)
    den <- sqrt(sum(a^2)) * sqrt(sum(b^2))
    if (den == 0) {
        return(NA_real_)
    }
    num / den
}

cosine_alignment_scores <- function(Gamma_hat, Gamma_true) {
    if (!is.matrix(Gamma_hat)) Gamma_hat <- matrix(Gamma_hat, ncol = 1L)
    if (!is.matrix(Gamma_true)) Gamma_true <- matrix(Gamma_true, ncol = 1L)
    apply(Gamma_hat, 2, function(g_hat) {
        max(apply(Gamma_true, 2, function(g_true) abs(cosine_similarity_r(g_hat, g_true))))
    })
}
