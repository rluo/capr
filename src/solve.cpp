#include <RcppArmadillo.h>

// ─────────────────────────────────────────────────────────────────────────────
//  β–STEP  (one Newton / Fisher like scoring update)
// ─────────────────────────────────────────────────────────────────────────────
//[[Rcpp::export]]
arma::vec newton_beta(const arma::cube& S, const arma::mat& X,
                      const arma::vec& T, const arma::vec& beta_init,
                      const arma::vec& gamma, int max_iter = 1,
                      double tol = 1e-6) {
  const arma::uword n = S.n_slices;
  arma::vec beta = beta_init;
  for (int it = 0; it < max_iter; ++it) {
    arma::mat Hbeta(beta.n_elem, beta.n_elem, arma::fill::zeros);
    arma::vec g(beta.n_elem, arma::fill::zeros);
    for (arma::uword i = 0; i < n; ++i) {
      double eta = arma::dot(X.row(i), beta);
      double wi = arma::trunc_exp(-eta) *
                  arma::as_scalar(gamma.t() * S.slice(i) * gamma);
      arma::vec xi = X.row(i).t();
      Hbeta += wi * (xi * xi.t()) / static_cast<double>(n);
      g += (T[i] - wi) * xi / static_cast<double>(n);
    }

#ifdef DEBUG
    std::cout << "Hbeta = \n" << Hbeta << std::endl;
#endif
    arma::vec delta = arma::solve(
        Hbeta, g, arma::solve_opts::fast + arma::solve_opts::likely_sympd);
#ifdef DEBUG
    std::cout << "delta = \n" << delta << std::endl;
#endif
    beta -= delta;  // step size 1 / ||delta||_inf to improve stability
    if (arma::norm(delta, "inf") < tol) {
      break;
    }
  }

#ifdef DEBUG
  std::cout << "beta = \n" << beta << std::endl;
#endif
  if (arma::norm(beta, "inf") > 16) {
    Rcpp::warning("Large beta values detected; results may be unstable.");
  }
  return beta;
}
