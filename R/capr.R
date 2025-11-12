#' Covariate Assisted Projection (CAP) Regression
#'
#' Fits (multiple) CAP components sequentially by alternating updates of
#' direction vectors \eqn{\gamma^{(k)}} and regression coefficients
#' \eqn{\beta^{(k)}}. Each component is estimated via a flip–flop algorithm,
#' with optional orthogonalization of successive directions.
#'
#' @param S Numeric 3D array of size \eqn{p \times p \times n} (stack of covariance matrices).
#' @param X Numeric matrix \eqn{n \times q} (design matrix).
#' @param weight Numeric vector of length \eqn{n} (default rep(1, n)).
#' @param Gamma.init Initial value of principal direction matrix  \eqn{\Gamma \in  R^{p\times k}} (default random Gaussian matrix).
#' @param B.init Initial value of coefficient \eqn{B \in R^{q \times K}} (default zero matrix).
#' @param K Integer scalar, number of components (\eqn{K \ge 1}).
#' @param max_iter Integer scalar, max flip–flop iterations per component (default 200).
#' @param tol Positive numeric scalar, convergence tolerance (default 1e-6).
#' @param orth Logical scalar; if TRUE (default), enforce orthogonality of successive \eqn{\gamma}.
#'
#' @return A list with:
#' \item{B}{numeric matrix \eqn{q \times K}, column \eqn{k} stores \eqn{\beta^{(k)}}}
#' \item{Gamma}{numeric matrix \eqn{p \times K}, column \eqn{k} stores \eqn{\gamma^{(k)}}}
#'
#' @details
#' For \eqn{k = 1, \dots, K}:
#' * Initialize \eqn{\beta^{(k)}} (zeros) and \eqn{\gamma^{(k)}} (random normal).
#' * If \code{orth = TRUE} and \code{k > 1}, orthogonalize \eqn{\gamma^{(k)}} against previous directions.
#' * Optionally rank-complete the covariance cube \code{S} given previous directions/coefficients
#'   if a function \code{rank_complete_s()} is available.
#' * Call \code{CAP_one_component()} to update \eqn{\beta^{(k)}} and \eqn{\gamma^{(k)}} until convergence.
#'
#' @note Requires \code{CAP_one_component()} to be available. If \code{rank_complete_s()} is not
#' defined, the input \code{S} is used as-is at each step.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' p <- 5
#' n <- 20
#' q <- 3
#' K <- 2
#' S <- array(0, dim = c(p, p, n))
#' for (i in 1:n) {
#'     A <- matrix(rnorm(p * p), p, p)
#'     S[, , i] <- crossprod(A) + diag(p)
#' }
#' X <- matrix(rnorm(n * q), n, q)
#' T <- rnorm(n)
#' # define CAP_one_component() in your package; then:
#' res <- capr(S, X, T, K)
#' str(res)
#' }
#' @seealso \code{\link{CAP_one_component}}, \code{\link{rank_complete_s}}
#' @export
capr <- function(S, X, K, Gamma.init = NULL, B.init = NULL, weight = NULL, max_iter = 200L, tol = 1e-6, orth = TRUE) {
    # ---- type checks ----
    # S: numeric p x p x n
    if (!is.array(S) || length(dim(S)) != 3L || !is.numeric(S)) {
        stop("`S` must be a numeric 3D array of shape p x p x n.", call. = FALSE)
    }

    dS <- dim(S)
    p <- dS[1]
    p2 <- dS[2]
    n <- dS[3]
    if (p != p2) stop("`S` must be square in its first two dims: dim(S)[1] == dim(S)[2].", call. = FALSE)

    # X: numeric n x q
    if (!is.matrix(X) || !is.numeric(X)) {
        stop("`X` must be a numeric matrix of shape n x q.", call. = FALSE)
    }
    if (nrow(X) != n) {
        stop(sprintf("nrow(X) must equal dim(S)[3] (= n = %d). Got nrow(X) = %d.", n, nrow(X)), call. = FALSE)
    }
    q <- ncol(X)

    # Scalars
    if (length(K) != 1L || !is.numeric(K) || is.na(K) || K < 1 || K != as.integer(K)) {
        stop("`K` must be a single positive integer (>= 1).", call. = FALSE)
    }
    K <- as.integer(K)

    if (is.null(Gamma.init)) {
        Gamma.init <- matrix(stats::rnorm(p * K), nrow = p, ncol = K)
    } else if (!is.matrix(Gamma.init) || !identical(dim(Gamma.init), c(p, K)) || !is.numeric(Gamma.init)) {
        stop("`Gamma.init` must be NULL, a numeric matrix of size p x K, or a numeric vector of length p*K.", call. = FALSE)
    }


    if (is.null(B.init)) {
        B.init <- matrix(0, nrow = q, ncol = K)
    } else if (!is.matrix(B.init) || !identical(dim(B.init), c(p, K)) || !is.numeric(B.init)) {
        stop("`B.init` must be NULL, a numeric matrix of size p x K, or a numeric vector of length p*K.", call. = FALSE)
    }


    if (is.null(weight)) weight <- rep(1, n)

    # weight: numeric vector length n
    if (!is.numeric(weight) || is.matrix(weight) || length(weight) != n) {
        stop(sprintf("`weight` must be a numeric vector of length n (= %d).", n), call. = FALSE)
    }


    if (length(max_iter) != 1L || !is.numeric(max_iter) || is.na(max_iter) ||
        max_iter < 1 || max_iter != as.integer(max_iter)) {
        stop("`max_iter` must be a single positive integer.", call. = FALSE)
    }
    max_iter <- as.integer(max_iter)

    if (length(tol) != 1L || !is.numeric(tol) || is.na(tol) || tol <= 0) {
        stop("`tol` must be a single positive number (> 0).", call. = FALSE)
    }

    if (length(orth) != 1L || !is.logical(orth) || is.na(orth)) {
        stop("`orth` must be a single logical (TRUE/FALSE).", call. = FALSE)
    }


    if (K == 1) {
        CAPre <- CAP_one_component(
            S, X, weight,
            beta_init = beta_k,
            gamma_init = gamma_k,
            Gprev = Gprev_step,
            max_iter = max_iter,
            tol = tol
        )
    } else {





    }







    # Optional routine
    has_rank_complete <- is.function(rank_complete_s)

    # ---- allocate outputs ----
    Gamma <- matrix(0, nrow = p, ncol = K) # columns are gamma^(k)
    B <- matrix(0, nrow = q, ncol = K) # columns are beta^(k)

    # Working copy of S (may be rank-completed each step)
    S_work <- S

    # Helper: orthogonalize v against columns of Gprev via QR of Gprev
    orthogonalise_qr <- function(v, Gprev) {
        stopifnot(is.numeric(v), length(v) == p)
        if (is.null(Gprev) || ncol(Gprev) == 0L) {
            return(as.numeric(v))
        }
        # QR of Gprev gives orthonormal basis Q for its column space
        qrobj <- qr(Gprev)
        Q <- tryCatch(qr.Q(qrobj), error = function(e) NULL)
        if (!is.null(Q)) {
            proj <- drop(Q %*% (crossprod(Q, v)))
            v <- v - proj
        } else {
            # Fallback to normal equations projection
            coef <- tryCatch(solve(crossprod(Gprev), crossprod(Gprev, v)),
                error = function(e) rep(0, ncol(Gprev))
            )
            v <- v - as.numeric(Gprev %*% coef)
        }
        # If nearly zero, jitter a bit
        nv <- sqrt(sum(v^2))
        if (!is.finite(nv) || nv < .Machine$double.eps) {
            v <- stats::rnorm(p)
            # remove Gprev component again
            if (!is.null(Q)) v <- v - drop(Q %*% (crossprod(Q, v)))
        }
        v
    }

    # ---- main loop over components ----
    for (k in seq_len(K)) {
        # initialize beta_k and gamma_k
        beta_k <- rep(0, q)
        gamma_k <- stats::rnorm(p)

        # If we already have previous directions, optionally orthogonalize and/or rank-complete
        if (k > 1L) {
            Gprev <- Gamma[, 1:(k - 1L), drop = FALSE]
            if (orth) {
                gamma_k <- orthogonalise_qr(gamma_k, Gprev)
            }
            if (has_rank_complete) {
                # rank_complete_s(S, X, Gprev, Bprev)
                Bprev <- B[, 1:(k - 1L), drop = FALSE]
                S_work <- rank_complete_s(S, X, Gprev, Bprev)
                # Basic sanity: keep dimensions
                if (!is.array(S_work) || length(dim(S_work)) != 3L ||
                    any(dim(S_work) != dim(S))) {
                    stop("`rank_complete_s()` must return a numeric array with the same dim as `S`.", call. = FALSE)
                }
            } else {
                S_work <- S
            }
        } else {
            S_work <- S
        }

        # Prepare Gprev_step to pass: either the matrix of previous gammas, or a 0-column matrix
        if (orth && k > 1L) {
            Gprev_step <- Gamma[, 1:(k - 1L), drop = FALSE]
        } else {
            Gprev_step <- matrix(numeric(0), nrow = p, ncol = 0L)
        }

        # One flip–flop update via user-supplied routine
        CAPre <- CAP_one_component(
            S_work, X, T,
            beta_init = beta_k,
            gamma_init = gamma_k,
            Gprev = Gprev_step,
            max_iter = max_iter,
            tol = tol
        )

        # Validate CAPre
        if (!is.list(CAPre) || !all(c("beta", "gamma") %in% names(CAPre))) {
            stop("`CAP_one_component()` must return a list with elements `beta` and `gamma`.", call. = FALSE)
        }
        beta_out <- CAPre$beta
        gamma_out <- CAPre$gamma
        if (!is.numeric(beta_out) || length(beta_out) != q) stop("Returned `beta` has wrong length.", call. = FALSE)
        if (!is.numeric(gamma_out) || length(gamma_out) != p) stop("Returned `gamma` has wrong length.", call. = FALSE)

        # store results
        B[, k] <- as.numeric(beta_out)
        Gamma[, k] <- as.numeric(gamma_out)
    }

    list(B = B, Gamma = Gamma)
}
