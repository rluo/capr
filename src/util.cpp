#include <RcppArmadillo.h>

#include "util.h"

//' Cosine similarity between two vectors
//'
//' Computes the cosine of the angle between \code{a} and \code{b}. Both
//' vectors must have the same length and non-zero norms.
//'
//' @noRd
// [[Rcpp::export]]
double cosine_similarity(const arma::vec& a, const arma::vec& b, double eps) {
  if (a.n_elem != b.n_elem) Rcpp::stop("cosine_similarity: size mismatch");

  const double na = arma::norm(a, 2);
  const double nb = arma::norm(b, 2);

  if (na < eps || nb < eps) Rcpp::stop("cosine_similarity: zero vector");

  return arma::dot(a, b) / (na * nb);
}
