#' Flury-Gautschi Common Principal Components
#'
#' Implements the Flury & Gautschi (1986) (FG) iterative algorithm and a variant to estimate a
#' common loading matrix across multiple covariance matrices. Each iteration
#' cycles over all ordered pairs of variable indices and updates a (2 x 2)
#' rotation so that the transformed matrices share diagonal structure.
#'
#' Two solvers are exported:
#' \describe{
#'   \item{`FG()`}{The original FG algorithm.}
#'   \item{`FG2()`}{An alternative algorithm by Eslami et al. (2013).}
#' }
#'
#' @param cov_array Numeric 3D array of shape \eqn{p x p x m}
#'   containing covariance matrices in its \eqn{m} slices.
#' @param p Optional integer specifying the matrix dimension; defaults to
#'   \code{dim(cov_array)[1]}.
#' @param m Optional integer specifying the number of matrices/slices; defaults
#'   to \code{dim(cov_array)[3]}.
#' @param maxit Integer scalar; number of outer iterations of the algorithm.
#'
#' @return A \eqn{p x p} numeric matrix of estimated common loadings.
#'
#' @references
#' Flury, B. N. (1984). Common principal components in k groups.
#' \emph{Journal of the American Statistical Association}, 79, 892-898.
#'
#' Flury, B. N., & Gautschi, W. (1986). An algorithm for simultaneous
#' orthogonal transformation of several positive definite symmetric matrices
#' to nearly diagonal form. \emph{SIAM Journal on Scientific and Statistical
#' Computing}, 7(1), 169-184.
#'
#' Eslami, A., Qannari, E. M., Kohler, A., & Bougeard, S. (2013).
#' General overview of methods of analysis of multi-group datasets.
#' \emph{Revue des Nouvelles Technologies de l'Information}, 25, 108-123.
#' @examples
#' \dontrun{
#' set.seed(1)
#' p <- 3
#' m <- 4
#' mats <- replicate(m,
#'     {
#'         A <- matrix(rnorm(p * p), p, p)
#'         crossprod(A)
#'     },
#'     simplify = FALSE
#' )
#' cov_cube <- array(NA_real_, dim = c(p, p, m))
#' for (k in 1:m) cov_cube[, , k] <- mats[[k]]
#' FG(cov_cube, maxit = 5)
#' FG2(cov_cube, maxit = 5)
#' }
#'
#' @export
FG <- function(cov_array, p = NULL, m = NULL, maxit = 30L) {
    args <- fg_validate_inputs(cov_array, p, m, maxit)
    weights <- rep(1, args$m)
    FG_cpp(args$cov_array, weights, args$maxit)
}

#' @rdname FG
#' @export
FG2 <- function(cov_array, p = NULL, m = NULL, maxit = 30L) {
    args <- fg_validate_inputs(cov_array, p, m, maxit)
    FG2_cpp(args$cov_array, args$maxit, args$p, args$m)
}

fg_validate_inputs <- function(cov_array, p, m, maxit) {
    if (!is.array(cov_array) || length(dim(cov_array)) != 3L || !is.numeric(cov_array)) {
        stop("`cov_array` must be a numeric p x p x m array.", call. = FALSE)
    }
    dims <- dim(cov_array)
    if (dims[1L] != dims[2L]) {
        stop("First two dimensions of `cov_array` must be equal (square matrices).", call. = FALSE)
    }
    if (is.null(p)) p <- dims[1L]
    if (is.null(m)) m <- dims[3L]

    if (length(maxit) != 1L || !is.numeric(maxit) || is.na(maxit) || maxit < 1) {
        stop("`maxit` must be a positive integer.", call. = FALSE)
    }
    if (length(p) != 1L || !is.numeric(p) || is.na(p) || p < 1) {
        stop("`p` must be a positive integer.", call. = FALSE)
    }
    if (length(m) != 1L || !is.numeric(m) || is.na(m) || m < 1) {
        stop("`m` must be a positive integer.", call. = FALSE)
    }

    if (dims[1L] != p || dims[3L] != m) {
        stop("Dimensions of `cov_array` must equal p x p x m.", call. = FALSE)
    }

    storage.mode(cov_array) <- "double"
    maxit <- as.integer(maxit)
    p <- as.integer(p)
    m <- as.integer(m)

    list(cov_array = cov_array, p = p, m = m, maxit = maxit)
}
