#' Covariate Assisted Principal (CAP) Regression
#'
#' Fits CAP components sequentially for principal direction vectors \eqn{\gamma^{(k)}} and regression
#' coefficients \eqn{\beta^{(k)}}, \eqn{k = 1, \ldots, K}. Each component is estimated via a flip-flop
#' algorithm with optional orthogonalization of successive directions.
#'
#' @param S Numeric 3D array of size \eqn{p \times p \times n} (for example, a stack of covariance matrices).
#' @param X Numeric matrix \eqn{n \times q} (design matrix), created for example by \code{model.matrix()}.
#' @param weight Numeric vector of length \eqn{n} (default \code{rep(1, n)}); each element should be proportional
#'   to the sample size for the corresponding slice \eqn{S[, , i]}.
#' @param Gamma.init Initial value of the principal direction array \eqn{\Gamma \in \mathbb{R}^{p \times n.init \times K}}
#'   (default: random Gaussian 3D array).
#' @param B.init Initial value of the coefficient array \eqn{B \in \mathbb{R}^{q \times n.init \times K}}
#'   (default: zero 3D array).
#' @param K Integer scalar; number of components (\eqn{K \ge 1}).
#' @param max_iter Integer scalar; maximum flip-flop iterations per component (default 200).
#' @param tol Positive numeric scalar; convergence tolerance (default 1e-6).
#' @param orth Logical scalar; if \code{TRUE} (default), enforce orthogonality of successive \eqn{\gamma^{(k)}}.
#'   If \code{FALSE}, no orthogonalization is performed (which may yield identical components).
#' @param n.init Integer scalar; number of random initializations (default 10). If \code{B.init} and
#'   \code{Gamma.init} are both supplied, \code{n.init} is ignored.
#' @return A list of class \code{capr} with:
#' \item{B}{numeric matrix \eqn{q \times K} whose \eqn{k}-th column stores \eqn{\beta^{(k)}}}
#' \item{Gamma}{numeric matrix \eqn{p \times K} whose \eqn{k}-th column stores \eqn{\gamma^{(k)}}}
#' \item{loglike}{negative log-likelihood, up to constant scaling and shift}
#' \item{S}{3D array used for fitting}
#' \item{X}{design matrix used for fitting}
#' \item{weight}{weight values used for fitting}
#'
#' @details
#' For component \eqn{k}, CAP solves
#' \deqn{\min_{\beta^{(k)}, \gamma^{(k)}} \frac{1}{2} \sum_{i=1}^{n} (\mathbf{x}_{i}^\top \beta^{(k)}) T_{i}
#' + \frac{1}{2} \sum_{i=1}^{n} \gamma^{(k)\top} S_{i} \gamma^{(k)} \exp(-\mathbf{x}_{i}^\top \beta^{(k)})}
#' subject to \deqn{\gamma^{(k)\top} H \gamma^{(k)} = 1} and, for \eqn{k > 1}, \deqn{\Gamma^{(k-1)\top} \gamma^{(k)} = \mathbf{0}.}
#' Here \eqn{T_{i}} denotes the weight for slice \eqn{i}, \eqn{S_{i}} is the \eqn{i}-th covariance slice,
#' and \eqn{H} is the positive definite matrix used for the orthogonality constraint (see Zhao et al., 2021).
#' The algorithm fits \eqn{\gamma^{(k)}} and \eqn{\beta^{(k)}} sequentially with multiple random initializations and
#' returns the solution pair that minimizes the negative log-likelihood.
#'
#' @references
#' Zhao, Y., Wang, B., Mostofsky, S. H., Caffo, B. S., & Luo, X. (2021).
#' "Covariate assisted principal regression for covariance matrix outcomes."
#' \emph{Biostatistics}, 22(3), 629-645.
#' @examples
#' simu.data <- simu.capr(seed = 123L, n = 120L)
#' K <- 2L
#' fit <- capr(
#'     S = simu.data$S,
#'     X = simu.data$X,
#'     K = K
#' )
#' print(fit)
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

    ## check if B.init and Gamma.init are properly specified
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
            B.init <- array(stats::rnorm(q * K * n.init), dim = c(q, n.init, K))
            Gamma.init <- array(stats::rnorm(p * K * n.init), dim = c(p, n.init, K))
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

    ret <- list(B = B_hat, Gamma = Gamma_hat, loglike = loglikevec, S = S, X = X, weight = weight)
    class(ret) <- c("capr")
    ret
}
