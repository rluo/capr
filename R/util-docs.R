#' Cosine similarity between numeric vectors
#'
#' Computes the cosine of the angle between two numeric vectors. Both vectors
#' must have equal length and non-zero Euclidean norms.
#'
#' @param a,b Numeric vectors of equal length.
#' @param eps Non-negative numeric tolerance used to guard against division by
#'   zero. Defaults to `1e-12`.
#'
#' @return A scalar double value in `[-1, 1]` representing the cosine
#'   similarity between `a` and `b`.
#' @examples
#' cosine_similarity(c(1, 2, 3), c(1, 2, 3))
#' cosine_similarity(c(1, 0), c(0, 1))
#' cosine_similarity(c(1, 2), c(-1, -2))
#' @export
cosine_similarity <- function(a, b, eps = 1e-12) {
    .Call(`_capr_cosine_similarity_cpp`, a, b, eps)
}

#' Log deviation from diagonality
#'
#' Evaluates the Flury-Gautschi log-deviation criterion for a collection of
#' covariance matrices transformed by a loading matrix.
#'
#' @param S_cube Numeric 3D array of shape `p x p x n` containing covariance
#'   matrices in its slices.
#' @param nval Numeric vector of length `n` giving weights for each matrix.
#' @param B Numeric `p x p` orthonormal matrix applied to the covariance slices.
#'
#' @return Numeric scalar value equal to
#'   \eqn{\sum_i n_i (\log \det \operatorname{diag}(B^\top S_i B)
#'   - \log \det(B^\top S_i B)) / (\sum_i n_i)}.
#' @examples
#' covs <- array(diag(2), dim = c(2, 2, 1))
#' log_deviation_from_diagonality(covs, 1, diag(2))
#' @export
log_deviation_from_diagonality <- function(S_cube, nval, B) {
    .Call(`_capr_log_deviation_from_diagonality_cpp`, S_cube, nval, B)
}
