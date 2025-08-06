# Load Rcpp to compile and load the C++ function
library(Rcpp)

# Compile the C++ file. This needs to be done once per R session.
# Make sure the .cpp file is in your working directory or provide the full path.
sourceCpp("dev/cap_g.cpp")

#' @title Covariate Assisted Principal (CAP) Regression - Algorithm 1
#'
#' @description
#' Implements the block coordinate descent algorithm from Zhao et al. (2021) to
#' compute the first pair of projection direction (gamma) and regression
#' coefficients (beta).
#'
#' @details
#' This function serves as a user-friendly wrapper for the core C++ implementation.
#' It automatically computes the sample covariance matrices `S_i`, the constraint
#' matrix `H` (as the average sample covariance), and provides random initial
#' values if none are supplied. The core computation is an iterative process
#' that alternates between updating `beta` via a Newton-Raphson step and
#' updating `gamma` by solving a generalized eigenvalue problem.
#'
#' @param Y A list of data matrices. Each element `Y[[i]]` is a `T_i x p` data
#'   matrix for subject `i`.
#' @param X An `n x q` matrix of covariate variables, where `n` is the number of
#'   subjects. It should include an intercept column if desired.
#' @param beta_init (Optional) A numeric vector of length `q` for the initial
#'   beta coefficients. If `NULL`, defaults to a vector of zeros.
#' @param gamma_init (Optional) A numeric vector of length `p` for the initial
#'   projection direction. If `NULL`, defaults to a random unit vector.
#' @param max_iter An integer specifying the maximum number of iterations.
#' @param tol A numeric value for the convergence tolerance. The algorithm stops
#'   when the Frobenius norm of the change in both beta and gamma is less than `tol`.
#'
#' @return A list with the following components:
#'   \item{beta}{The estimated `q x 1` regression coefficient vector.}
#'   \item{gamma}{The estimated `p x 1` projection direction vector.}
#'   \item{iterations}{The number of iterations performed before convergence or stopping.}
#'   \item{objective}{The final value of the objective function (negative log-likelihood, ignoring constants) to be minimized.}
#'
#' @author Gemini
#' @references
#' Zhao, Y., Wang, B., Mostofsky, S. H., Caffo, B. S., & Luo, X. (2021).
#' Covariate Assisted Principal regression for covariance matrix outcomes.
#' Biostatistics, 22(3), 629-645.
#'
#' @export
#' @examples
#' # --- Generate Sample Data ---
#' # Set parameters
#' n <- 50 # number of subjects
#' p <- 10 # number of variables
#' q <- 2 # number of covariates (e.g., intercept and one variable)
#' Ti <- 30 # number of observations per subject
#'
#' set.seed(123)
#'
#' # True parameters
#' true_gamma <- rnorm(p)
#' true_gamma <- true_gamma / sqrt(sum(true_gamma^2))
#' true_beta <- c(0.5, -1.5)
#'
#' # Generate data
#' X <- cbind(1, rnorm(n))
#' Y_list <- vector("list", n)
#'
#' for (i in 1:n) {
#'     log_var <- X[i, ] %*% true_beta
#'     variance <- exp(log_var)
#'     # Simplified simulation: scale a base covariance
#'     cov_matrix <- diag(p) + tcrossprod(true_gamma) * as.numeric(variance)
#'     # Ensure positive definiteness
#'     if (min(eigen(cov_matrix)$values) <= 0) {
#'         cov_matrix <- cov_matrix + diag(p) * 1e-4
#'     }
#'     Y_list[[i]] <- MASS::mvrnorm(n = Ti, mu = rep(0, p), Sigma = cov_matrix)
#' }
#'
#' # --- Run the Algorithm ---
#' result <- cap_algorithm1(Y = Y_list, X = X)
#'
#' # --- View Results ---
#' print("Estimated Beta:")
#' print(result$beta)
#'
#' print("True Beta:")
#' print(true_beta)
#'
#' # Note: Estimated gamma may be -true_gamma (sign is arbitrary)
#' print("Inner product of estimated and true gamma (should be close to 1 or -1):")
#' print(crossprod(result$gamma, true_gamma))
#'
#' print(paste("Converged in", result$iterations, "iterations."))
cap_algorithm1 <- function(Y, X, beta_init = NULL, gamma_init = NULL, max_iter = 100, tol = 1e-6) {
    # ---- 1. Input Validation ----
    if (!is.list(Y)) stop("Y must be a list of matrices.")
    if (!is.matrix(X)) stop("X must be a matrix.")

    n <- length(Y)
    p <- ncol(Y[[1]])
    q <- ncol(X)

    if (n != nrow(X)) stop("Number of subjects in Y (length of list) must match number of rows in X.")
    if (any(sapply(Y, ncol) != p)) stop("All matrices in Y must have the same number of columns (p).")

    # ---- 2. Parameter Initialization ----
    if (is.null(beta_init)) {
        beta_init <- rep(0, q)
    }

    if (is.null(gamma_init)) {
        gamma_init <- rnorm(p)
        gamma_init <- gamma_init / sqrt(sum(gamma_init^2))
    }

    # ---- 3. Pre-computation of H ----
    # As suggested in the paper (Section 3), H is the average sample covariance matrix.
    # H = sum(S_i) / sum(T_i), where S_i = Y_i' * Y_i

    S_list <- lapply(Y, function(y) crossprod(y))
    T_vec <- sapply(Y, nrow)

    S_sum <- Reduce("+", S_list)
    T_sum <- sum(T_vec)

    H_matrix <- S_sum / T_sum

    # ---- 4. Call the C++ Backend ----

    result <- cap_algorithm1_cpp(
        Y_list = Y,
        X = X,
        H = H_matrix,
        beta_init = as.vector(beta_init),
        gamma_init = as.vector(gamma_init),
        max_iter = as.integer(max_iter),
        tol = as.double(tol)
    )

    # ---- 5. Calculate Final Objective Value and Format Output ----
    # Objective function L = 1/2 * sum(T_i * (x_i' * beta)) + 1/2 * sum(exp(-x_i' * beta) * gamma' * S_i * gamma)

    final_beta <- result$beta
    final_gamma <- result$gamma

    term1 <- 0.5 * sum(T_vec * (X %*% final_beta))

    term2_sum <- 0
    exp_vals <- exp(-X %*% final_beta)
    for (i in 1:n) {
        gamma_s_gamma <- crossprod(final_gamma, S_list[[i]]) %*% final_gamma
        term2_sum <- term2_sum + exp_vals[i] * gamma_s_gamma
    }

    result$objective <- term1 + 0.5 * term2_sum

    # Set names for nicer output
    names(result$beta) <- colnames(X)
    names(result$gamma) <- colnames(Y[[1]])

    return(result)
}
