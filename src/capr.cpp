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
  double loglike;
};

// ─────────────────────────────────────────────────────────────────────────────
//   One CAP component – flip-flop,    returning  CAPResult
// ─────────────────────────────────────────────────────────────────────────────
static CAPResult CAP_one_component_core(const arma::cube& S, const arma::mat& X,
                                        const arma::vec& T,
                                        const arma::mat& beta_init,
                                        const arma::mat& gamma_init,
                                        OptGamma opt_Gamma_prev,
                                        int max_iter = 200, double tol = 1e-6) {
  //  matrix  with columns β and γ, for m different initializations
  int max_inits = beta_init.n_cols;

  arma::mat gamma_work = gamma_init;  // copy, so we can normalise columns

  arma::vec best_gamma = gamma_init.col(0) * 0;
  arma::vec best_beta = beta_init.col(0) * 0;
  double best_loglike = arma::datum::inf;

  const arma::uword p = S.n_rows, n = S.n_slices;

  //  H  =  average Σ S_i
  arma::mat H(p, p, arma::fill::zeros);
  for (arma::uword i = 0; i < n; ++i) H += S.slice(i);
  H /= static_cast<double>(n);

  for (arma::uword i = 0; i < gamma_work.n_cols; ++i) {
    arma::vec gamma = gamma_work.col(i);

    // normalize to H  norm = 1
    // --- current quadratic form  γᵀ H γ  -----------------------------------
    const double qf = arma::as_scalar(gamma.t() * H * gamma);

    if (qf <= tol)
      Rcpp::stop("gamma has (near-)zero H-norm — cannot normalise");

    // --- scale γ so that  γᵀ H γ = 1  ---------------------------------------
    gamma /= std::sqrt(qf);

    gamma_work.col(i) = gamma;
  }

  const bool has_constraints =
      opt_Gamma_prev && opt_Gamma_prev->get().n_cols > 0;

  for (int itinit = 0; itinit < max_inits; ++itinit) {
    arma::vec beta = beta_init.col(itinit);
    arma::vec gamma = gamma_work.col(itinit);

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

    double loglike = cap_loglike_cpp(S, X, T, beta, gamma);

    if (loglike < best_loglike) {
      best_loglike = loglike;
      best_beta = beta;
      best_gamma = gamma;
    }
  }
  return {std::move(best_beta), std::move(best_gamma), best_loglike};  // NRVO
}

// [[Rcpp::export]]
Rcpp::List CAP_one_component_unconstrained(
    const arma::cube& S, const arma::mat& X, const arma::vec& T,
    const arma::mat beta0, const arma::mat gamma0, int max_iter = 200,
    double tol = 1e-6) {
  CAPResult res = CAP_one_component_core(S, X, T, beta0, gamma0, std::nullopt,
                                         max_iter, tol);

  return Rcpp::List::create(Rcpp::Named("beta") = res.beta,
                            Rcpp::Named("gamma") = res.gamma,
                            Rcpp::Named("loglike") = res.loglike);
}

// [[Rcpp::export]]
Rcpp::List CAP_one_component(const arma::cube& S, const arma::mat& X,
                             const arma::vec& T, const arma::mat beta0,
                             const arma::mat gamma0,
                             const arma::mat& Gamma_prev, int max_iter = 200,
                             double tol = 1e-6) {
  CAPResult res = CAP_one_component_core(S, X, T, beta0, gamma0,
                                         std::cref(Gamma_prev), max_iter, tol);

  return Rcpp::List::create(Rcpp::Named("beta") = res.beta,
                            Rcpp::Named("gamma") = res.gamma,
                            Rcpp::Named("loglike") = res.loglike);
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
    const arma::cube& Binit, const arma::cube& Gammainit,
    const bool orth = true, const int max_iter = 200, const double tol = 1e-6) {
  // Binit and Gammainit are  variable x  initializations x components
  const arma::uword p = S.n_rows;
  const arma::uword q = X.n_cols;

  arma::mat Gamma(p, K, arma::fill::zeros);  // p × K
  arma::mat B(q, K, arma::fill::zeros);      // q × K
  arma::vec loglikevec(K);                   // K × 1
  loglikevec.fill(arma::datum::inf);
  for (int k = 0; k < K; ++k) {
    arma::mat beta_k = Binit.slice(k);       // q x m
    arma::mat gamma_k = Gammainit.slice(k);  // p x m

    arma::cube S_current = S;
    OptGamma gamma_prev_opt = std::nullopt;
    std::optional<arma::mat> Gprev_holder;

    if (k > 0) {
      Gprev_holder.emplace(Gamma.cols(0, k - 1));
      const arma::mat& Gprev = *Gprev_holder;
      arma::mat Bprev = B.cols(0, k - 1);

      S_current = rank_complete_s(S, X, Gprev, Bprev);

      if (orth) {
        for (arma::uword init_idx = 0; init_idx < gamma_k.n_cols; ++init_idx) {
          gamma_k.col(init_idx) =
              orthogonalise_qr(gamma_k.col(init_idx), Gprev);
        }
        gamma_prev_opt = std::cref(Gprev);
      }
    }

    CAPResult res = CAP_one_component_core(S_current, X, T, beta_k, gamma_k,
                                           gamma_prev_opt, max_iter, tol);

    Gamma.col(k) = res.gamma;
    B.col(k) = res.beta;
    loglikevec(k) = res.loglike;
  }

  return Rcpp::List::create(Rcpp::Named("B") = B, Rcpp::Named("Gamma") = Gamma,
                            Rcpp::Named("loglike") = loglikevec);
}
