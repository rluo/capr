#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

/*
 * cov_cube: P x P x M covariance array
 * maxit: number of iterations
 * P: number of variables (dimension)
 * M: number of subjects/covariance matrices
 * returns: P x P matrix of common loadings
 */

// [[Rcpp::export]]
arma::mat FG_cpp(const arma::cube& cov_cube, int maxit, int P, int M) {
  const arma::uword p_dim = static_cast<arma::uword>(P);
  const arma::uword m_dim = static_cast<arma::uword>(M);

  if (cov_cube.n_rows != p_dim || cov_cube.n_cols != p_dim ||
      cov_cube.n_slices != m_dim) {
    Rcpp::stop("FG_cpp: `cov_cube` must be a P x P x M array.");
  }

  arma::mat B = arma::eye<arma::mat>(p_dim, p_dim);

  for (int it = 0; it < maxit; ++it) {
    for (int p = 0; p < P - 1; ++p) {
      for (int e = p + 1; e < P; ++e) {
        arma::mat Q = arma::eye<arma::mat>(2, 2);
        arma::mat Mmat = arma::eye<arma::mat>(2, 2);
        arma::uvec idx = {static_cast<arma::uword>(p),
                          static_cast<arma::uword>(e)};
        arma::mat H = B.cols(idx);

        for (arma::uword k = 0; k < m_dim; ++k) {
          arma::mat Tk = H.t() * cov_cube.slice(k) * H;
          double d1 = arma::as_scalar(Q.col(0).t() * Tk * Q.col(0));
          double d2 = arma::as_scalar(Q.col(1).t() * Tk * Q.col(1));
          Mmat += ((d1 - d2) / (d1 * d2)) * Tk;
        }

        arma::vec eigval;
        arma::mat eigvec;
        arma::eig_sym(eigval, eigvec, Mmat);

        Q = eigvec;
        H = B.cols(idx);
        B.cols(idx) = H * Q;
      }
    }
  }
  return B;
}
