#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' @title Core CAP Algorithm 1 in C++
//' @description
//' Implements the block coordinate descent algorithm described in Algorithm 1
//' of Zhao et al. (2021) to find the first pair of projection direction (gamma)
//' and regression coefficients (beta).
//'
//' @param Y_list A list of data matrices, where each element is a T_i x p
// matrix for subject i. ' @param X A numeric matrix of covariates (n x q). '
//@param H A p x p positive definite matrix for the constraint on gamma. '
//@param beta_init Initial values for the beta coefficients (q x 1 vector). '
//@param gamma_init Initial values for the gamma projection vector (p x 1
// vector). ' @param max_iter Maximum number of iterations. ' @param tol
// Convergence tolerance. ' @return A list containing the estimated 'beta',
//'gamma', and the number of 'iterations'. ' @keywords internal
// [[Rcpp::export]]
Rcpp::List cap_algorithm1_cpp(const Rcpp::List& Y_list, const arma::mat& X,
                              const arma::mat& H, arma::vec beta_init,
                              arma::vec gamma_init, int max_iter, double tol) {
  // ---- 1. Data Preparation ----

  int n = Y_list.size();  // Number of subjects
  int p = H.n_rows;       // Number of variables
  int q = X.n_cols;       // Number of covariates

  // Store S_i = Y_i' * Y_i and T_i = #rows(Y_i)
  arma::field<arma::mat> S_field(n);
  arma::vec T_vec(n);

  for (int i = 0; i < n; ++i) {
    arma::mat Y_i = Rcpp::as<arma::mat>(Y_list[i]);
    S_field(i) = Y_i.t() * Y_i;
    T_vec(i) = static_cast<double>(Y_i.n_rows);
  }

  arma::vec beta = beta_init;
  arma::vec gamma = gamma_init;
  int iter;

  // ---- 2. Main Iteration Loop ----

  for (iter = 0; iter < max_iter; ++iter) {
    arma::vec beta_old = beta;
    arma::vec gamma_old = gamma;

    // ---- Step (1): Update beta given gamma ----
    // This is a Newton-Raphson step: beta_new = beta - H^-1 * g
    // where g is the gradient and H is the Hessian of the objective function.

    arma::vec gradient = arma::zeros<arma::vec>(q);
    arma::mat hessian = arma::zeros<arma::mat>(q, q);

    for (int i = 0; i < n; ++i) {
      arma::vec xi = X.row(i).t();

      // Scalar term: exp(-x_i' * beta) * (gamma' * S_i * gamma)
      double gamma_s_gamma = arma::as_scalar(gamma.t() * S_field(i) * gamma);
      double exp_val = std::exp(arma::as_scalar(-X.row(i) * beta));

      // Update gradient
      gradient += (T_vec(i) - exp_val * gamma_s_gamma) * xi;

      // Update Hessian
      hessian += (exp_val * gamma_s_gamma) * (xi * xi.t());
    }

    // Solve the linear system and update beta
    beta =
        beta - arma::solve(hessian, gradient, arma::solve_opts::likely_sympd);

    // ---- Step (2): Update gamma given beta ----
    // Minimize (1/2) * gamma' * A * gamma subject to gamma' * H * gamma = 1
    // where A = sum_i( exp(-x_i' * beta) * S_i )
    // This is a generalized eigenvalue problem: A * v = lambda * H * v
    // The solution is the eigenvector 'v' corresponding to the smallest
    // eigenvalue 'lambda'.

    arma::mat A = arma::zeros<arma::mat>(p, p);
    for (int i = 0; i < n; ++i) {
      double exp_val = std::exp(arma::as_scalar(-X.row(i) * beta));
      A += exp_val * S_field(i);
    }

    arma::vec eigval;
    arma::mat eigvec;

    // eig_pair finds generalized eigenvalues and sorts them in ascending order.
    // The solution for gamma is the eigenvector for the smallest eigenvalue.
    bool success = arma::eig_pair(eigval, eigvec, A, H);

    if (!success) {
      Rcpp::stop("Generalized eigenvalue decomposition failed.");
    }

    gamma = eigvec.col(0);  // Smallest is the first column

    // ---- 3. Convergence Check ----

    double beta_change = arma::norm(beta - beta_old, "fro");
    // Handle potential sign flipping in eigenvector
    double gamma_change = std::min(arma::norm(gamma - gamma_old, "fro"),
                                   arma::norm(gamma + gamma_old, "fro"));

    if (beta_change < tol && gamma_change < tol) {
      iter++;  // Increment to reflect the final iteration count
      break;
    }
  }

  if (iter == max_iter) {
    Rcpp::warning("Algorithm did not converge after %d iterations.", max_iter);
  }

  // ---- 4. Return Results ----

  return Rcpp::List::create(Rcpp::Named("beta") = beta,
                            Rcpp::Named("gamma") = gamma,
                            Rcpp::Named("iterations") = iter);
}