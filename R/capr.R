#' Covariate Assisted Projection (CAP) Regression
#'
#' Fits (multiple) CAP components sequentially by alternating updates of
#' direction vectors \eqn{\gamma^{(k)}} and regression coefficients
#' \eqn{\beta^{(k)}}. Each component is estimated via a flip–flop algorithm,
#' with optional orthogonalization of successive directions.
#'
#' @param S Numeric 3D array of size \eqn{p \times p \times n} (stack of covariance matrices).
#' @param X Numeric matrix \eqn{n \times q} (design matrix).
#' @param weight Numeric vector of length \eqn{n} (default rep(1, n)), where each element should be proportional the sample size for each \eqn{S_i}.
#' @param Gamma.init Initial value of principal direction matrix  \eqn{\Gamma \in  R^{p\times n.init \times K}} (default random Gaussian cube or array).
#' @param B.init Initial value of coefficient \eqn{B \in R^{q \times n.init \times K}} (default zero cube or array).
#' @param K Integer scalar, number of components (\eqn{K \ge 1}).
#' @param max_iter Integer scalar, max flip–flop iterations per component (default 200).
#' @param tol Positive numeric scalar, convergence tolerance (default 1e-6).
#' @param orth Logical scalar; if TRUE (default), enforce orthogonality of successive \eqn{\gamma}.
#' @param n.init Integer scalar; number of random initializations (default 10). If B.init and Gamma.init are set properly, n.init is ignored.
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
#' simu.data <- simu.capr(seed = 123L, p = 5L, n = 120L)
#' K <- 2L
#' fit <- capr(
#'     S = simu.data$S,
#'     X = simu.data$X,
#'     K = K
#' )
#' print(fit)
#' @seealso \code{\link{CAP_one_component}}, \code{\link{rank_complete_s}}
#' @export
capr <- function(S, X, K, B.init = NULL, Gamma.init = NULL, weight = NULL, max_iter = 200L, tol = 1e-6, orth = TRUE, n.init = 10L) {
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

    ## check  if B.init and Gamma.init are properly specifieed
    if (!is.null(B.init) && !is.null(Gamma.init)) {
        if (is.numeric(B.init) && is.numeric(Gamma.init) && dim(B.init)[3L] == dim(Gamma.init)[3L] &&
            dim(B.init)[2L] == dim(Gamma.init)[2L]) {
            n.init <- dim(B.init)[2L]
            if (length(B.init) != q * n.init * K || length(Gamma.init) != p * n.init * K) {
                stop("If both `B.init` and `Gamma.init` are specified, they must be numeric arrays of shape (q x n.init x K) and (p x n.init x K), respectively.", call. = FALSE)
            }
            B.init <- array(B.init, dim = c(q, n.init, K))
            Gamma.init <- array(Gamma.init, dim = c(p, n.init, K))
        } else {
            stop("If both `B.init` and `Gamma.init` are specified, they must be numeric arrays of shape (q x n.init x K) and (p x n.init x K), respectively.", call. = FALSE)
        }
    } else {
        if (!is.null(n.init)) {
            B.init <- array(rnorm(q * K * n.init), dim = c(q, n.init, K))
            Gamma.init <- array(rnorm(p * K * n.init), dim = c(p, n.init, K))
        } else {
            stop("If `B.init` and `Gamma.init` are not both specified, `n.init` must be provided.", call. = FALSE)
        }
    }

    storage.mode(B.init) <- "double"
    storage.mode(Gamma.init) <- "double"
    if (is.null(weight)) weight <- rep(1, n)
    # weight: numeric vector length n
    if (!is.numeric(weight) || is.matrix(weight) || length(weight) != n) {
        stop(sprintf("`weight` must be a numeric vector of length n (= %d).", n), call. = FALSE)
    }
    weight <- as.numeric(weight)
    if (any(!is.finite(weight))) {
        stop("`weight` must contain only finite values.", call. = FALSE)
    }

    if (any(weight <= 0)) {
        stop("`weight` must contain only positive values.", call. = FALSE)
    }

    weight <- weight / max(weight)
    storage.mode(weight) <- "double"

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

    if (!orth && K > 1L) {
        warning("Fitting multiple components without orthogonalization may lead to identical components.", call. = FALSE)
    }

    cap_fit <- CAP_multi_components(
        S = S,
        X = X,
        T = weight,
        K = K,
        Binit = B.init,
        Gammainit = Gamma.init,
        orth = orth,
        max_iter = max_iter,
        tol = tol
    )

    B_hat <- cap_fit$B
    Gamma_hat <- cap_fit$Gamma
    loglikevec <- cap_fit$loglike

    if (!is.null(colnames(X))) {
        rownames(B_hat) <- colnames(X)
    }
    colnames(B_hat) <- paste0("Comp", seq_len(K))
    colnames(Gamma_hat) <- paste0("Comp", seq_len(K))
    rownames(Gamma_hat) <- paste0("V", seq_len(p))

    ret <- list(B = B_hat, Gamma = Gamma_hat, loglike = loglikevec)
    class(ret) <- c("capr")
    ret
}
