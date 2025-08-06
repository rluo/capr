// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <cmath>
using namespace Rcpp;
using namespace arma;

// Function to compute sample covariance matrices S_i
// [[Rcpp::export]]
arma::cube compute_sample_covariances(const List& Y) {
  int n = Y.size();
  arma::mat Y_first = as<arma::mat>(Y[0]);
  int p = Y_first.n_cols;

  arma::cube S(p, p, n);

  for (int i = 0; i < n; i++) {
    arma::mat Y_i = as<arma::mat>(Y[i]);
    S.slice(i) = Y_i.t() * Y_i;
  }

  return S;
}

// Function to compute H matrix (average sample covariance)
// [[Rcpp::export]]
arma::mat compute_H_matrix(const arma::cube& S, const arma::vec& T_vec) {
  int n = S.n_slices;
  int p = S.n_rows;

  arma::mat H = arma::zeros<arma::mat>(p, p);
  double total_T = arma::sum(T_vec);

  for (int i = 0; i < n; i++) {
    H += S.slice(i) / total_T;
  }

  return H;
}

// Function to solve generalized eigenvalue problem for gamma update
// [[Rcpp::export]]
arma::vec solve_gamma_update(const arma::cube& S, const arma::mat& X,
                             const arma::vec& beta, const arma::mat& H) {
  int n = S.n_slices;
  int p = S.n_rows;

  // Compute weighted sum of covariance matrices
  arma::mat A = arma::zeros<arma::mat>(p, p);

  for (int i = 0; i < n; i++) {
    double weight = std::exp(-arma::dot(X.row(i), beta));
    A += weight * S.slice(i);
  }

  // Solve generalized eigenvalue problem: A * gamma = lambda * H * gamma
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_gen(eigval, eigvec, A, H);

  // Find the eigenvector corresponding to the smallest eigenvalue
  arma::uword min_idx;
  eigval.min(min_idx);

  arma::vec gamma = arma::real(eigvec.col(min_idx));

  // Normalize to satisfy constraint gamma^T H gamma = 1
  double norm_factor = std::sqrt(arma::dot(gamma, H * gamma));
  gamma = gamma / norm_factor;

  return gamma;
}

// Function to compute beta update using Newton-Raphson
// [[Rcpp::export]]
arma::vec compute_beta_update(const arma::cube& S, const arma::mat& X,
                              const arma::vec& T_vec, const arma::vec& beta,
                              const arma::vec& gamma) {
  int n = S.n_slices;
  int q = X.n_cols;

  arma::mat hessian = arma::zeros<arma::mat>(q, q);
  arma::vec gradient = arma::zeros<arma::vec>(q);

  for (int i = 0; i < n; i++) {
    arma::vec x_i = X.row(i).t();
    double exp_term = std::exp(-arma::dot(x_i, beta));
    double quad_form = arma::dot(gamma, S.slice(i) * gamma);

    // Compute gradient contribution
    gradient += (T_vec(i) - exp_term * quad_form) * x_i;

    // Compute Hessian contribution
    hessian += exp_term * quad_form * (x_i * x_i.t());
  }

  // Newton-Raphson update: beta_new = beta_old - H^(-1) * g
  arma::vec beta_update = beta - arma::solve(hessian, gradient);

  return beta_update;
}

// Main Algorithm 1 implementation
// [[Rcpp::export]]
List cap_algorithm_cpp(const List& Y, const arma::mat& X,
                       const arma::vec& T_vec, const arma::vec& beta_init,
                       const arma::vec& gamma_init, double tol = 1e-6,
                       int max_iter = 1000) {
  // Compute sample covariance matrices
  arma::cube S = compute_sample_covariances(Y);

  // Compute H matrix
  arma::mat H = compute_H_matrix(S, T_vec);

  // Initialize parameters
  arma::vec beta = beta_init;
  arma::vec gamma = gamma_init;

  // Normalize initial gamma to satisfy constraint
  double norm_factor = std::sqrt(arma::dot(gamma, H * gamma));
  gamma = gamma / norm_factor;

  // Main iteration loop
  for (int iter = 0; iter < max_iter; iter++) {
    arma::vec beta_old = beta;
    arma::vec gamma_old = gamma;

    // Update beta
    beta = compute_beta_update(S, X, T_vec, beta, gamma);

    // Update gamma
    gamma = solve_gamma_update(S, X, beta, H);

    // Check convergence
    double beta_diff = arma::norm(beta - beta_old);
    double gamma_diff = arma::norm(gamma - gamma_old);

    if (beta_diff < tol && gamma_diff < tol) {
      Rcout << "Converged after " << iter + 1 << " iterations" << std::endl;
      break;
    }

    if (iter == max_iter - 1) {
      Rcout << "Maximum iterations reached without convergence" << std::endl;
    }
  }

  // Compute final objective value
  double obj_val = 0.0;
  for (int i = 0; i < S.n_slices; i++) {
    arma::vec x_i = X.row(i).t();
    double quad_form = arma::dot(gamma, S.slice(i) * gamma);
    obj_val += 0.5 * arma::dot(x_i, beta) * T_vec(i) +
               0.5 * quad_form * std::exp(-arma::dot(x_i, beta));
  }

  return List::create(Named("beta") = beta, Named("gamma") = gamma,
                      Named("objective") = obj_val, Named("S") = S,
                      Named("H") = H);
}
