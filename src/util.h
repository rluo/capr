#pragma once
#include <RcppArmadillo.h>

/**
 * @title Cosine Similarity Between Two Vectors
 * @name cosine_similarity
 * @description
 * Computes the cosine similarity between two numeric vectors
 * \eqn{a} and \eqn{b}:
 * \deqn{\cos(\theta) \;=\; \frac{a^\top b}{\|a\|_2 \; \|b\|_2}}
 * where \eqn{\|\cdot\|_2} denotes the Euclidean (L2) norm.
 *
 * @details
 * - The function accepts any Armadillo dense vector type (e.g. `arma::vec`,
 *   `arma::rowvec`, or subviews).
 * - Both vectors must be the same length; otherwise an error is thrown.
 * - If either vector has (near-)zero norm (smaller than \code{eps}),
 *   the function throws an error.
 *
 * @param a A numeric vector.
 * @param b A numeric vector of the same length as \code{a}.
 * @param eps A small numeric tolerance (default \code{1e-12}) used to guard
 *   against division by zero if one of the vectors is near zero length.
 *
 * @return A scalar numeric value in \eqn{[-1, 1]} giving the cosine similarity.
 *
 * @examples
 * # Two identical vectors => cosine similarity = 1
 * cosine_similarity(c(1, 2, 3), c(1, 2, 3))
 *
 * # Orthogonal vectors => cosine similarity = 0
 * cosine_similarity(c(1, 0), c(0, 1))
 *
 * # Opposite vectors => cosine similarity = -1
 * cosine_similarity(c(1, 2), c(-1, -2))
 *
 * @seealso \code{\link[=Rcpp]{Rcpp}}, \code{\link{RcppArmadillo}}
 * @export
 */
// [[Rcpp::export]]
static inline double cosine_similarity(const arma::vec& a, const arma::vec& b,
                                       const double eps = 1e-12) {
  if (a.n_elem != b.n_elem) Rcpp::stop("cosine_similarity: size mismatch");

  const double na = arma::norm(a, 2);
  const double nb = arma::norm(b, 2);

  if (na < eps || nb < eps) Rcpp::stop("cosine_similarity: zero vector");

  return arma::dot(a, b) / (na * nb);
}