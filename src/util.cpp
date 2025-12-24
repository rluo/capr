#include "util.h"

#include <RcppArmadillo.h>

//' Cosine similarity between two vectors
//'
//' Computes the cosine of the angle between \code{a} and \code{b}. Both
//' vectors must have the same length and non-zero norms.
//'
//' @noRd
// [[Rcpp::export(name = "cosine_similarity_cpp")]]
double cosine_similarity_cpp(const arma::vec& a, const arma::vec& b,
                             double eps) {
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
// [[Rcpp::export(name = "log_deviation_from_diagonality_cpp")]]
double log_deviation_from_diagonality_cpp(const arma::cube& S_cube,
                                          const arma::vec& nval,
                                          const arma::mat& B) {
  if (S_cube.n_slices != nval.n_elem)
    Rcpp::stop("S_cube third dimension (%d) != length(nval) (%d)",
               S_cube.n_slices, nval.n_elem);
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
  sum = sum / arma::accu(nval);
  return sum;
}

// Compute the (twice) CAP log-likelihood
// \frac{1}{2}\sum_{i=1}^{n}\Bigl(\mathbf{x}_{i}^{\top}\boldsymbol{\beta}\Bigr)\,T_{i}\;+\;\frac{1}{2}\sum_{i=1}^{n}\boldsymbol{\gamma}^{\top}\mathbf{S}_{i}^{(k)}\boldsymbol{\gamma}\;\exp\!\bigl(-\mathbf{x}_{i}^{\top}\boldsymbol{\beta}\bigr)
// [[Rcpp::export(name = "cap_loglike_cpp")]]
double cap_loglike_cpp(const arma::cube& S_cube, const arma::mat& X,
                       const arma::vec& T, const arma::vec& beta,
                       const arma::vec& gamma) {
  size_t n = S_cube.n_slices;
  double loglike = 0.0;
  arma::vec XB = X * beta;  // n x 1
  for (size_t i = 0; i < n; ++i) {
    double BTAB = arma::as_scalar(gamma.t() * S_cube.slice(i) * gamma);
    loglike += XB[i] * T[i] + BTAB * std::exp(-XB[i]);
  }
  return loglike;
}

// [[Rcpp::export]]
arma::mat GtSG(const arma::mat& G, const arma::cube& S) {
  // const arma::uword p = G.n_rows;
  const arma::uword K = G.n_cols;
  const arma::uword n = S.n_slices;

  arma::mat result(n, K, arma::fill::zeros);

  for (arma::uword i = 0; i < K; ++i) {
    arma::vec g = G.col(i);
    for (arma::uword j = 0; j < n; ++j) {
      result(j, i) = arma::as_scalar(g.t() * S.slice(j) * g);
    }
  }
  return result;
}