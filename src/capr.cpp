// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <algorithm>
#include <functional>
#include <optional>

#include "complete.h"
#include "solve.h"
#include "util.h"

using OptGamma = std::optional<std::reference_wrapper<const arma::mat>>;

struct CAPResult {
  arma::vec beta;
  arma::vec gamma;
};

// ─────────────────────────────────────────────────────────────────────────────
//   One CAP component – flip-flop,  now returning  CAPResult
// ─────────────────────────────────────────────────────────────────────────────

static CAPResult CAP_one_component_core(const arma::cube& S, const arma::mat& X,
                                        const arma::vec& T,
                                        const arma::vec& beta_init,
                                        const arma::vec& gamma_init,
                                        OptGamma opt_Gamma_prev,
                                        int max_iter = 200, double tol = 1e-6) {
  arma::vec beta = beta_init;
  arma::vec gamma = gamma_init;

  const arma::uword p = S.n_rows, n = S.n_slices;

  //  H  =  average Σ S_i
  arma::mat H(p, p, arma::fill::zeros);
  for (arma::uword i = 0; i < n; ++i) H += S.slice(i);
  H /= static_cast<double>(n);

  // normalize to H  norm = 1
  // --- current quadratic form  γᵀ H γ  -----------------------------------
  const double qf = arma::as_scalar(gamma.t() * H * gamma);

  if (qf <= tol) Rcpp::stop("gamma has (near-)zero H-norm — cannot normalise");

  // --- scale γ so that  γᵀ H γ = 1  ---------------------------------------
  gamma /= std::sqrt(qf);

  const bool has_constraints =
      opt_Gamma_prev && opt_Gamma_prev->get().n_cols > 0;

  for (int it = 0; it < max_iter; ++it) {
    arma::vec beta_old = beta, gamma_old = gamma;

    beta = newton_beta(S, X, T, beta, gamma);

    // build A(β)
    arma::mat A(p, p, arma::fill::zeros);
    for (arma::uword i = 0; i < n; ++i)
      A += std::exp(-arma::dot(X.row(i), beta)) * S.slice(i);

    gamma = has_constraints ? solve_gamma(A, H, opt_Gamma_prev->get())
                            : solve_gamma_unconstrained(A, H);

    if (std::max(arma::norm(beta - beta_old, "inf"),
                 arma::norm(gamma - gamma_old, "inf")) < tol)
      break;
  }

  return {std::move(beta), std::move(gamma)};  // NRVO
}

// [[Rcpp::export]]
Rcpp::List CAP_one_component_unconstrained(const arma::cube& S,
                                           const arma::mat& X,
                                           const arma::vec& T, arma::vec beta0,
                                           arma::vec gamma0, int max_iter = 200,
                                           double tol = 1e-6) {
  CAPResult res = CAP_one_component_core(S, X, T, beta0, gamma0, std::nullopt,
                                         max_iter, tol);

  return Rcpp::List::create(Rcpp::Named("beta") = res.beta,
                            Rcpp::Named("gamma") = res.gamma);
}

// [[Rcpp::export]]
Rcpp::List CAP_one_component(const arma::cube& S, const arma::mat& X,
                             const arma::vec& T, arma::vec beta0,
                             arma::vec gamma0, const arma::mat& Gamma_prev,
                             int max_iter = 200, double tol = 1e-6) {
  CAPResult res = CAP_one_component_core(S, X, T, beta0, gamma0,
                                         std::cref(Gamma_prev), max_iter, tol);

  return Rcpp::List::create(Rcpp::Named("beta") = res.beta,
                            Rcpp::Named("gamma") = res.gamma);
}

arma::vec orthogonalise_qr(const arma::vec& gamma_k,
                           const arma::mat& Gamma_prev) {
  if (Gamma_prev.n_cols == 0) return normalise(gamma_k);

  arma::mat Q_R;
  arma::mat Rtrash;
  arma::qr_econ(Q_R, Rtrash, Gamma_prev);  // Q_R: p×m, orthonormal

  arma::vec proj = Q_R.t() * gamma_k;  // m×1
  arma::vec g = normalise(gamma_k - Q_R * proj);
  return g;
}

// [[Rcpp::export]]
Rcpp::List CAP_multi_components(
    const arma::cube& S, const arma::mat& X, const arma::vec& T, const int K,
    const arma::mat& Binit, const arma::mat& Gammainit, const bool orth = true,
    const int max_iter = 200, const double tol = 1e-6) {
  const arma::uword p = S.n_rows;
  const arma::uword q = X.n_cols;

  arma::mat Gamma(p, K, arma::fill::zeros);
  arma::mat B(q, K, arma::fill::zeros);

  for (int k = 0; k < K; ++k) {
    arma::vec beta_k = Binit.col(k);
    arma::vec gamma_k = Gammainit.col(k);

    arma::cube S_current = S;
    OptGamma gamma_prev_opt = std::nullopt;
    std::optional<arma::mat> Gprev_holder;

    if (k > 0) {
      Gprev_holder.emplace(Gamma.cols(0, k - 1));
      const arma::mat& Gprev = *Gprev_holder;
      arma::mat Bprev = B.cols(0, k - 1);

      S_current = rank_complete_s(S, X, Gprev, Bprev);

      if (orth) {
        gamma_k = orthogonalise_qr(gamma_k, Gprev);
        gamma_prev_opt = std::cref(Gprev);
      }
    }

    CAPResult res = CAP_one_component_core(S_current, X, T, beta_k, gamma_k,
                                           gamma_prev_opt, max_iter, tol);

    Gamma.col(k) = res.gamma;
    B.col(k) = res.beta;
  }

  return Rcpp::List::create(Rcpp::Named("B") = B, Rcpp::Named("Gamma") = Gamma);
}
