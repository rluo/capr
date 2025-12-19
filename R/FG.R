#' Flury-Gautschi Common Principal Components
#'
#' Implements the Flury & Gautschi (1986) iterative algorithm to estimate a
#' common loading matrix across multiple covariance matrices. Each iteration
#' cycles over all ordered pairs of variable indices and updates a \(2 x 2\)
#' rotation so that the transformed matrices share diagonal structure.
#'
#' Two solvers are exported:
#' \describe{
#'   \item{`FG()`}{Calls the original implementation in \code{FG.cpp}.}
#'   \item{`FG2()`}{Calls the alternative implementation in \code{FG2.cpp} that
#'   follows the same interface but can be used to benchmark or validate the
#'   computation.}
#' }
#'
#' @param cov_array Numeric 3D array of shape \eqn{P x P x M}
#'   containing covariance matrices in its \eqn{M} slices.
#' @param P Optional integer specifying the matrix dimension; defaults to
#'   \code{dim(cov_array)[1]}.
#' @param M Optional integer specifying the number of matrices/slices; defaults
#'   to \code{dim(cov_array)[3]}.
#' @param maxit Integer scalar; number of outer iterations of the algorithm.
#'
#' @return A \eqn{P x P} numeric matrix of estimated common loadings.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' P <- 3
#' M <- 4
#' mats <- replicate(M,
#'     {
#'         A <- matrix(rnorm(P * P), P, P)
#'         crossprod(A)
#'     },
#'     simplify = FALSE
#' )
#' cov_cube <- array(NA_real_, dim = c(P, P, M))
#' for (k in 1:M) cov_cube[, , k] <- mats[[k]]
#' FG(cov_cube, maxit = 5)
#' FG2(cov_cube, maxit = 5)
#' }
#'
#' @seealso \code{\link{CAP_one_component}}, \code{\link{rank_complete_s}}
#' @export
FG <- function(cov_array, P = NULL, M = NULL, maxit = 30L) {
    args <- fg_validate_inputs(cov_array, P, M, maxit)
    weights <- rep(1, args$M)
    FG_cpp(args$cov_array, weights, args$maxit)
}

#' @rdname FG
#' @export
FG2 <- function(cov_array, P = NULL, M = NULL, maxit = 30L) {
    args <- fg_validate_inputs(cov_array, P, M, maxit)
    FG2_cpp(args$cov_array, args$maxit, args$P, args$M)
}

fg_validate_inputs <- function(cov_array, P, M, maxit) {
    if (!is.array(cov_array) || length(dim(cov_array)) != 3L || !is.numeric(cov_array)) {
        stop("`cov_array` must be a numeric P x P x M array.", call. = FALSE)
    }
    dims <- dim(cov_array)
    if (dims[1L] != dims[2L]) {
        stop("First two dimensions of `cov_array` must be equal (square matrices).", call. = FALSE)
    }
    if (is.null(P)) P <- dims[1L]
    if (is.null(M)) M <- dims[3L]

    if (length(maxit) != 1L || !is.numeric(maxit) || is.na(maxit) || maxit < 1) {
        stop("`maxit` must be a positive integer.", call. = FALSE)
    }
    if (length(P) != 1L || !is.numeric(P) || is.na(P) || P < 1) {
        stop("`P` must be a positive integer.", call. = FALSE)
    }
    if (length(M) != 1L || !is.numeric(M) || is.na(M) || M < 1) {
        stop("`M` must be a positive integer.", call. = FALSE)
    }

    if (dims[1L] != P || dims[3L] != M) {
        stop("Dimensions of `cov_array` must equal P x P x M.", call. = FALSE)
    }

    storage.mode(cov_array) <- "double"
    maxit <- as.integer(maxit)
    P <- as.integer(P)
    M <- as.integer(M)

    list(cov_array = cov_array, P = P, M = M, maxit = maxit)
}
