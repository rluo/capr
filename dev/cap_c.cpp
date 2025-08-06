// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// -----------------------------------------------------------------------------
//  Helper: build S_i = Y_iᵀ Y_i and record T_i  (row-count of Y_i)
// -----------------------------------------------------------------------------
static void make_covariances(const Rcpp::List& Y_list,
                             arma::field<arma::mat>& S_list, arma::vec& T_vec) {
  const std::size_t n = Y_list.size();
  S_list.set_size(n);
  T_vec.set_size(n);

  for (std::size_t i = 0; i < n; ++i) {
    arma::mat Yi = Rcpp::as<arma::mat>(Y_list[i]);  // T_i × p
    T_vec(i) = Yi.n_rows;                           // T_i
    S_list(i) = Yi.t() * Yi;                        // p × p
  }
}

// -----------------------------------------------------------------------------
//  Core optimiser implementing Algorithm 1
// -----------------------------------------------------------------------------
Rcpp::List cap_first_direction_core(const arma::field<arma::mat>& S_list,
                                    const arma::vec& T_vec,
                                    const arma::mat& X,  // n × q
                                    arma::vec beta,      // q
                                    arma::vec gamma,     // p
                                    const double tol = 1e-6,
                                    const unsigned max_iter = 1000) {
  const std::size_t n = X.n_rows;
  const std::size_t q = X.n_cols;
  const std::size_t p = gamma.n_rows;

  // ---------------------------------------------------------------------
  // H  :=  pooled sample covariance  (eq. 3.1, paper)
  // ---------------------------------------------------------------------
  arma::mat H(p, p, arma::fill::zeros);
  double total_T = arma::accu(T_vec);

  for (std::size_t i = 0; i < n; ++i) H += S_list(i);
  H /= total_T;

  // ---------------------------------------------------------------------
  // Iterative BCD
  // ---------------------------------------------------------------------
  double obj_old = arma::datum::inf, obj_new = 0.0, diff = tol + 1.0;
  unsigned iter = 0;

  while (iter < max_iter && diff > tol) {
    //------------------------------------------------------------------
    // (1)   β-update  (Newton step – eq. 3.2)
    //------------------------------------------------------------------
    arma::mat H_beta(q, q, arma::fill::zeros);
    arma::vec g_beta(q, arma::fill::zeros);

    for (std::size_t i = 0; i < n; ++i) {
      const double expo = std::exp(-arma::dot(X.row(i), beta));
      const double gSg = arma::as_scalar(gamma.t() * S_list(i) * gamma);
      const double w = expo * gSg;  // weight
      const arma::rowvec xi = X.row(i);

      H_beta += w * xi.t() * xi;          // Hessian
      g_beta += (T_vec(i) - w) * xi.t();  // gradient
    }
    beta += arma::solve(H_beta, g_beta, arma::solve_opts::fast);

    //------------------------------------------------------------------
    // (2)   γ-update  (Prop. A.1 ⇒ smallest gen. eigen-pair)
    //------------------------------------------------------------------
    arma::mat A(p, p, arma::fill::zeros);
    for (std::size_t i = 0; i < n; ++i) {
      const double expo = std::exp(-arma::dot(X.row(i), beta));
      A += expo * S_list(i);
    }
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A, H);  // generalised
    gamma = eigvec.col(0);                // minimise ⇒ λ_min

    // renormalise so that γᵀ H γ = 1  (constraint)
    gamma /= std::sqrt(arma::as_scalar(gamma.t() * H * gamma));

    //------------------------------------------------------------------
    // Objective value
    //------------------------------------------------------------------
    obj_new = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      const double lin = arma::dot(X.row(i), beta);  // x_iᵀ β
      const double gSg = arma::as_scalar(gamma.t() * S_list(i) * gamma);
      obj_new += 0.5 * (lin * T_vec(i) + gSg * std::exp(-lin));
    }

    diff = std::abs(obj_old - obj_new);
    obj_old = obj_new;
    ++iter;
  }

  return Rcpp::List::create(
      Rcpp::Named("beta") = beta, Rcpp::Named("gamma") = gamma,
      Rcpp::Named("objective") = obj_old, Rcpp::Named("iterations") = iter,
      Rcpp::Named("converged") = (diff <= tol));
}

// -----------------------------------------------------------------------------
//  [[Rcpp::export]]  –  User-facing wrapper
//      Y_list : list of Ti×p numeric matrices (one per subject)
//      X      : n×q numeric matrix (1st column should be the intercept)
// -----------------------------------------------------------------------------
[[Rcpp::export]]
Rcpp::List cap_first_direction_cpp(Rcpp::List Y_list, const arma::mat& X,
                                   arma::vec beta, arma::vec gamma,
                                   const double tol = 1e-6,
                                   const unsigned max_iter = 1000) {
  arma::field<arma::mat> S_list;
  arma::vec T_vec;
  make_covariances(Y_list, S_list, T_vec);

  return cap_first_direction_core(S_list, T_vec, X, beta, gamma, tol, max_iter);
}
