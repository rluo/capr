#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

/*
 * cov_Group: a list of M PxP covariance matrices
 * L: number of iterations
 * P: number of variables (dimension)
 * M: number of groups
 * returns: PxP matrix of common loadings
 */

// [[Rcpp::export]]
arma::mat FG_cpp(const Rcpp::List& cov_Group, int L, int P, int M) {
  arma::mat B = arma::eye<arma::mat>(P, P);

  for (int l = 0; l < L; ++l) {
    for (int p = 0; p < P - 1; ++p) {
      for (int e = p + 1; e < P; ++e) {
        arma::mat Q = arma::eye<arma::mat>(2, 2);
        arma::mat Mmat = arma::eye<arma::mat>(2, 2);
        std::vector<arma::mat> T(M);
        arma::vec d1(M);
        arma::vec d2(M);

        for (int k = 0; k < M; ++k) {
          arma::mat Ck = Rcpp::as<arma::mat>(cov_Group[k]);
          arma::mat H = B.cols(arma::uvec({(unsigned)p, (unsigned)e}));
          T[k] = H.t() * Ck * H;
          d1(k) = arma::as_scalar(Q.col(0).t() * T[k] * Q.col(0));
          d2(k) = arma::as_scalar(Q.col(1).t() * T[k] * Q.col(1));
          Mmat += (d1(k) - d2(k)) / (d1(k) * d2(k)) * T[k];
        }

        arma::vec eigval;
        arma::mat eigvec;
        arma::eig_sym(eigval, eigvec, Mmat);

        Q = eigvec;  // columns are eigenvectors
        arma::mat H = B.cols(arma::uvec({(unsigned)p, (unsigned)e}));
        B.cols(arma::uvec({(unsigned)p, (unsigned)e})) = H * Q;
      }
    }
  }
  return B;
}
