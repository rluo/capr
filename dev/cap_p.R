#' Covariate Assisted Principal (CAP) Regression Algorithm
#'
#' Implementation of Algorithm 1 from the CAP regression paper
#'
#' @param Y List of data matrices, where each element is a T_i x p matrix
#' @param X n x q matrix of covariate variables (first column should be ones for intercept)
#' @param beta_init Initial values for beta (q x 1 vector)
#' @param gamma_init Initial values for gamma (p x 1 vector)
#' @param tol Convergence tolerance (default: 1e-6)
#' @param max_iter Maximum number of iterations (default: 1000)
#' @return List containing estimated beta, gamma, objective value, and other components


library(RcppArmadillo)
library(Rcpp)
# Source the C++ functions
sourceCpp("src/cap_p.cpp")

#' Main CAP Algorithm wrapper function
cap_algorithm <- function(Y, X, beta_init = NULL, gamma_init = NULL,
                          tol = 1e-6, max_iter = 1000) {
    # Input validation
    if (!is.list(Y)) {
        stop("Y must be a list of matrices")
    }

    if (!is.matrix(X)) {
        stop("X must be a matrix")
    }

    n <- length(Y)
    if (nrow(X) != n) {
        stop("Number of rows in X must match length of Y")
    }

    # Get dimensions
    p <- ncol(Y[[1]])
    q <- ncol(X)

    # Check that all Y matrices have the same number of columns
    if (!all(sapply(Y, ncol) == p)) {
        stop("All matrices in Y must have the same number of columns")
    }

    # Compute T_vec (number of observations per subject)
    T_vec <- sapply(Y, nrow)

    # Initialize parameters if not provided
    if (is.null(beta_init)) {
        beta_init <- rep(0, q)
        beta_init[1] <- log(mean(T_vec)) # Initialize intercept
    }

    if (is.null(gamma_init)) {
        # Initialize gamma as the first eigenvector of the average covariance
        S_avg <- Reduce("+", lapply(1:n, function(i) t(Y[[i]]) %*% Y[[i]])) / sum(T_vec)
        eigen_decomp <- eigen(S_avg)
        gamma_init <- eigen_decomp$vectors[, 1]
    }

    # Call C++ implementation
    result <- cap_algorithm_cpp(Y, X, T_vec, beta_init, gamma_init, tol, max_iter)

    # Add additional R-specific components
    result$call <- match.call()
    result$n <- n
    result$p <- p
    result$q <- q
    result$T_vec <- T_vec

    class(result) <- "cap_fit"
    return(result)
}

#' Print method for cap_fit objects
print.cap_fit <- function(x, ...) {
    cat("CAP Regression Results\n")
    cat("=====================\n\n")
    cat("Dimensions:\n")
    cat(sprintf("  n (subjects): %d\n", x$n))
    cat(sprintf("  p (variables): %d\n", x$p))
    cat(sprintf("  q (covariates): %d\n", x$q))
    cat(sprintf("  Total observations: %d\n", sum(x$T_vec)))
    cat("\n")

    cat("Regression coefficients (beta):\n")
    print(x$beta)
    cat("\n")

    cat("Projection direction (gamma):\n")
    print(x$gamma)
    cat("\n")

    cat(sprintf("Final objective value: %.6f\n", x$objective))
}

#' Summary method for cap_fit objects
summary.cap_fit <- function(object, ...) {
    cat("CAP Regression Summary\n")
    cat("=====================\n\n")

    # Print basic info
    print(object)

    # Additional diagnostics
    cat("\nDiagnostics:\n")
    cat(sprintf(
        "  Gamma norm (should be ~1): %.6f\n",
        sqrt(sum(object$gamma^2))
    ))

    # Check constraint satisfaction: gamma^T H gamma = 1
    constraint_val <- t(object$gamma) %*% object$H %*% object$gamma
    cat(sprintf(
        "  Constraint satisfaction (gamma^T H gamma): %.6f\n",
        constraint_val
    ))
}

#' Generate synthetic data for testing
generate_cap_data <- function(n = 50, p = 5, q = 2, T_i = 100,
                              beta_true = c(1, -0.5), gamma_true = NULL) {
    # Generate gamma if not provided
    if (is.null(gamma_true)) {
        gamma_true <- rnorm(p)
        gamma_true <- gamma_true / sqrt(sum(gamma_true^2))
    }

    # Generate covariates
    X <- cbind(1, rnorm(n)) # Intercept and one covariate

    # Generate covariance matrices
    Sigma_list <- list()
    Y_list <- list()

    for (i in 1:n) {
        # Generate eigenvalues based on log-linear model
        log_lambda <- X[i, ] %*% beta_true
        lambda_i <- exp(log_lambda)

        # Create covariance matrix (simplified model)
        # In practice, this would be more complex
        D <- diag(c(lambda_i, rep(1, p - 1)))

        # Generate random orthogonal matrix
        Q <- qr.Q(qr(matrix(rnorm(p * p), p, p)))

        # Ensure gamma_true is (approximately) an eigenvector
        Q[, 1] <- gamma_true
        Q <- qr.Q(qr(Q)) # Re-orthogonalize

        Sigma_i <- Q %*% D %*% t(Q)
        Sigma_list[[i]] <- Sigma_i

        # Generate data
        Y_i <- mvrnorm(T_i, mu = rep(0, p), Sigma = Sigma_i)
        Y_list[[i]] <- Y_i
    }

    return(list(
        Y = Y_list,
        X = X,
        Sigma_list = Sigma_list,
        beta_true = beta_true,
        gamma_true = gamma_true
    ))
}

#' Example usage and testing
test_cap_algorithm <- function() {
    library(MASS) # for mvrnorm

    cat("Testing CAP Algorithm Implementation\n")
    cat("===================================\n\n")

    # Generate test data
    set.seed(123)
    data <- generate_cap_data(n = 30, p = 4, q = 2, T_i = 50)

    cat("Generated synthetic data:\n")
    cat(sprintf("  n = %d subjects\n", length(data$Y)))
    cat(sprintf("  p = %d variables\n", ncol(data$Y[[1]])))
    cat(sprintf("  q = %d covariates\n", ncol(data$X)))
    cat("\n")

    cat("True parameters:\n")
    cat("Beta (true):", data$beta_true, "\n")
    cat("Gamma (true):", round(data$gamma_true, 3), "\n\n")

    # Fit CAP model
    cat("Fitting CAP model...\n")
    fit <- cap_algorithm(data$Y, data$X, tol = 1e-6, max_iter = 500)

    # Print results
    cat("\nResults:\n")
    print(fit)

    # Compare with true values
    cat("\nComparison with true values:\n")
    cat("Beta estimation error:", round(fit$beta - data$beta_true, 4), "\n")
    cat("Gamma correlation with true:", round(abs(cor(fit$gamma, data$gamma_true)), 4), "\n")

    return(fit)
}
