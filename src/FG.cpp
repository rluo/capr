// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "util.h"

// Canonicalize: order by decreasing eigenvalue, flip sign so first entry ≥ 0
inline void canonicalize_Q(arma::mat& Q, const arma::vec& eigval) {
  // If eigval(0) < eigval(1), swap cols 0 and 1
  if (eigval(0) < eigval(1)) {
    Q.swap_cols(0, 1);
  }
  // For each column, flip sign so first nonzero entry is positive
  for (arma::uword j = 0; j < Q.n_cols; ++j) {
    arma::uword idx = arma::abs(Q.col(j)).index_max();  // new syntax
    if (Q(idx, j) < 0) {
      Q.col(j) *= -1;
    }
  }
}

// Computes the maximum columnwise norm difference between A and B, up to sign
// flip. Returns max_j min(norm(A_j - B_j), norm(A_j + B_j)) for all columns j
double max_col_diff_up_to_sign(const arma::mat& A, const arma::mat& B) {
  if (A.n_rows != B.n_rows || A.n_cols != B.n_cols)
    throw std::invalid_argument("Matrix size mismatch");

  double max_diff = 0.0;
  for (arma::uword j = 0; j < A.n_cols; ++j) {
    double d1 = arma::norm(A.col(j) - B.col(j), "inf");
    double d2 = arma::norm(A.col(j) + B.col(j), "inf");
    max_diff = std::max(max_diff, std::min(d1, d2));
  }
  return max_diff;
}

arma::mat FG(const arma::cube& arr, const arma::vec& nval, int max_iter = 100,
             double epsilon = 1e-6) {
  int M = arr.n_slices;
  int P = arr.n_rows;
  arma::mat B = arma::eye(P, P);
  // arma::mat B_old = B;
  double crit_old;
  // std::vector<arma::mat> T(M, arma::zeros(2, 2));
  arma::cube T(2, 2, M, arma::fill::zeros);
  // std::vector<double> d1(M, 0.0), d2(M, 0.0);
  arma::vec d1(M, arma::fill::zeros), d2(M, arma::fill::zeros);

  for (int l = 0; l < max_iter; ++l) {
    crit_old = log_deviation_from_diagonality_cpp(arr, nval, B);
    // B_old = B;
    for (int p = 0; p < P - 1; ++p) {
      for (int e = p + 1; e < P; ++e) {
        arma::mat Q = arma::eye(2, 2);
        // arma::mat Q = arma::mat(2, 2, arma::fill::randu);
        arma::mat Q_old = Q;
        arma::mat H(P, 2);  // allocate once
        arma::mat update(P, 2);

        for (int iter = 0; iter < max_iter; ++iter) {
          arma::mat Mmat = arma::zeros(2, 2);

          for (int k = 0; k < M; ++k) {
            const arma::mat& Ck = arr.slice(k);  // avoid copy
            H.col(0) = B.col(p);
            H.col(1) = B.col(e);
            T.slice(k) = H.t() * Ck * H;
            const arma::vec& Q0 = Q.col(0);
            const arma::vec& Q1 = Q.col(1);
            d1[k] = arma::as_scalar(Q0.t() * T.slice(k) * Q0);
            d2[k] = arma::as_scalar(Q1.t() * T.slice(k) * Q1);
            Mmat +=
                (nval[k] * ((d1[k] - d2[k]) / (d1[k] * d2[k]))) * T.slice(k);
          }
          arma::vec eigval;
          arma::mat eigvec;
          arma::eig_sym(eigval, eigvec, Mmat, "std");
          canonicalize_Q(eigvec, eigval);
          Q_old = Q;
          Q = eigvec;

          // Stopping criterion: change in Q below threshold
          if (max_col_diff_up_to_sign(Q, Q_old) < epsilon) {
// #ifndef NDEBUG
#if 0
            std::cout << "Q early return: " << iter << " " << l
                      << " p e = " << p << " " << e << "\n";
            std::cout << "Mmat = \n" << Mmat << "\n";

#endif
            break;
          }
        }
        // #ifndef NDEBUG
#if 0
        std::cout << l << ": Q diff = " << max_col_diff_up_to_sign(Q, Q_old)
                  << "\n";
        std::cout << "Q = \n" << Q << "\n";
        std::cout << "Q_old = \n" << Q_old << "\n";
#endif
        H.col(0) = B.col(p);
        H.col(1) = B.col(e);
        update = H * Q;
        B.col(p) = update.col(0);
        B.col(e) = update.col(1);
      }
    }
    if (crit_old - log_deviation_from_diagonality_cpp(arr, nval, B) < epsilon) {
      break;
    }
  }
  return B;
}

// [[Rcpp::export]]
arma::mat FG_cpp(const arma::cube& S_cube,  // p × p × n   (numeric array)
                 const arma::vec& nval,     // length n    (numeric vector)
                 const int max_iter = 100, const double epsilon = 1e-6) {
  if (S_cube.n_slices != nval.n_elem)
    Rcpp::stop("S_cube third dimension (%d) != length(nval) (%d)",
               S_cube.n_slices, nval.n_elem);

  return FG(S_cube, nval, max_iter, epsilon);  // ← call your routine
}
