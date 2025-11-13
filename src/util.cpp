#include "util.h"

#include <RcppArmadillo.h>

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

//' Log deviation from diagonality
//'
//' Computes \eqn{\sum_i n_i (\log \det \operatorname{diag}(B^\top S_i B) -
//' \log \det(B^\top S_i B))}, the criterion minimized by the FG algorithm.
//'
//' @noRd
// [[Rcpp::export]]
double log_deviation_from_diagonality(const arma::cube& S_cube,
                                      const arma::vec& nval,
                                      const arma::mat& B) {
  assert(S_cube.n_slices == nval.size());
  double sum = 0.0;
  for (size_t i = 0; i < S_cube.n_slices; ++i) {
    arma::mat BTAB = B.t() * S_cube.slice(i) * B;
    arma::vec diag_elems = BTAB.diag();
    // Check for non-positive values before log, to avoid nan
    if (arma::any(diag_elems <= 0)) {
      throw std::runtime_error("Non-positive diagonal element in BTAB.diag()");
    }
    double logdet_diag = arma::sum(arma::log(diag_elems));
    double sign = 0.0;
    double logdet_BTAB = 0.0;
    arma::log_det(logdet_BTAB, sign, BTAB);
    sum += nval[i] * (logdet_diag - logdet_BTAB);
  }
  return sum;
}
