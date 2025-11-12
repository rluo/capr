// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <functional>
#include <optional>

#include "complete.hpp"
#include "solve.hpp"
#include "util.hpp"

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

  const has_constraints = opt_Gamma_pre && opt_Gamma_pre->get().n_cols > 0;

  for (int it = 0; it < max_iter; ++it) {
    arma::vec beta_old = beta, gamma_old = gamma;

    beta = newton_beta(S, X, T, beta, gamma);

    // build A(β)
    arma::mat A(p, p, arma::fill::zeros);
    for (arma::uword i = 0; i < n; ++i)
      A += std::exp(-arma::dot(X.row(i), beta)) * S.slice(i);

    if (has_constraints) {
      gamma = solve_gamma(A, H, opt_Gamma_prev->get());
    } else {
      gamma = solve_gamma(A, H);
    }

    if (std::max(arma::norm(beta - beta_old, "inf"),
                 arma::norm(gamma - gamma_old, "inf")) < tol)
      break;
  }

  return {std::move(beta), std::move(gamma)};  // NRVO
}

// ─────────────────────────────────────────────────────────────────────────────
//  Rcpp export  (calls the new value-returning routine)
// ─────────────────────────────────────────────────────────────────────────────
//' One-component CAP update  (β–γ flip-flop)
//'
//' @param S          3-D array (p × p × n) with slices S_i
//' @param X          Design matrix (n × q)
//' @param T          Length-n vector
//' @param beta0      Numeric start vector
//' @param gamma0     Numeric start vector
//' @param Gamma_prev p × (k-1) matrix (may have 0 columns)
//' @param max_iter   integer, default 200
//' @param tol        double, default 1e-6
//' @return           list(beta = .., gamma = ..)
// [[Rcpp::export]]
Rcpp::List CAP_one_component(const arma::cube& S, const arma::mat& X,
                             const arma::vec& T, arma::vec beta0,
                             arma::vec gamma0, const arma::mat& Gamma_prev,
                             int max_iter = 200, double tol = 1e-6) {
  CAPResult res =
      CAP_one_component_core(S, X, T, beta0, gamma0, Gamma_prev, max_iter, tol);

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

/**
 * @title multi-Component Covariate Assisted Projection (CAP) Regression
 * @name capr
 * @description
 * Fits multiple CAP components sequentially by alternating updates of
 * direction vectors \eqn{\gamma^{(k)}} and regression coefficients
 * \eqn{\beta^{(k)}}. Each component is estimated using a flip–flop algorithm,
 * with optional orthogonalization of successive directions.
 * @details
 * For \eqn{k = 1, \dots, K}:
 * - Initialize \eqn{\beta^{(k)}} (zeros) and \eqn{\gamma^{(k)}} (random
 * normal).
 * - If \code{k > 0} and \code{orth = TRUE}, orthogonalize \eqn{\gamma^{(k)}}
 *   against previously found directions.
 * - Optionally rank-complete the covariance cube \eqn{S} given previous
 *   directions and coefficients.
 * - Run \code{CAP_one_component()} until convergence or until \code{max_iter}
 *   is reached, updating \eqn{\beta^{(k)}} and \eqn{\gamma^{(k)}}.
 * - Store results in the matrices \code{B} and \code{Gamma}.
 *
 * @param S A numeric cube of size \eqn{p \times p \times n}, representing
 *   a stack of covariance matrices \eqn{S_i}.
 * @param X A numeric matrix of size \eqn{n \times q}, the design matrix
 *   (rows \eqn{x_i^\top}).
 * @param T A numeric vector of length \eqn{n}, auxiliary data such as
 *   time indices or response variable depending on the CAP variant.
 * @param K Integer, number of CAP components to fit.
 * @param max_iter Maximum number of flip–flop iterations per component
 *   (default \code{200}).
 * @param tol Numeric convergence tolerance (default \code{1e-6}).
 * @param orth Logical; if \code{TRUE} (default), enforce orthogonality
 *   of successive direction vectors by QR-orthogonalization.
 *
 * @return An R list with two components:
 * \item{B}{Matrix of regression coefficients of size \eqn{q \times K}.
 *   Column \eqn{k} stores \eqn{\beta^{(k)}}.}
 * \item{Gamma}{Matrix of direction vectors of size \eqn{p \times K}.
 *   Column \eqn{k} stores \eqn{\gamma^{(k)}}.}
 *
 * @note
 * - Requires \code{CAP_one_component()} and \code{rank_complete_s()} to be
 *   available in the package namespace.
 * - Random initialization of \eqn{\gamma^{(k)}} may yield slightly
 *   different solutions across runs.
 *
 * @examples
 * \dontrun{
 *   set.seed(123)
 *   p <- 5; n <- 20; q <- 3; K <- 2
 *   S <- array(0, dim = c(p, p, n))
 *   for (i in 1:n) {
 *     A <- matrix(rnorm(p*p), p, p)
 *     S[,,i] <- crossprod(A) + diag(p)
 *   }
 *   X <- matrix(rnorm(n*q), n, q)
 *   T <- rnorm(n)
 *   res <- capr(S, X, T, K)
 *   str(res)
 * }
 *
 * @seealso \code{\link{CAP_one_component}}, \code{\link{rank_complete_s}}
 * @export
 */
// [[Rcpp::export]]
Rcpp::List capr(const arma::cube& S, const arma::mat& X, const arma::vec& T,
                const int K, const int max_iter = 200, const double tol = 1e-6,
                const bool orth = true) {
  const arma::uword p = S.n_rows;
  const arma::uword q = X.n_cols;

  arma::mat Gamma(p, K, arma::fill::zeros);  // store γ^{(1)},…,γ^{(K)}
  arma::mat B(q, K, arma::fill::zeros);      // store β^{(1)},…,β^{(K)}

  arma::mat Gamma_prev(p, K, arma::fill::zeros);  // pre-allocated

  arma::cube S_work = S;  // may be rank-completed each step

  for (int k = 0; k < K; ++k) {
    // initialize β and γ
    arma::vec beta_k(q, arma::fill::zeros);
    arma::vec gamma_k = arma::randn<arma::vec>(p);
    // gamma_k /= arma::norm(gamma_k, 2);

    // rank-complete S if orthogonality step required
    if (k > 0) {
      // use first k columns of Gamma_prev
      arma::mat Gprev = Gamma.cols(0, k - 1);
      gamma_k = orthogonalise_qr(gamma_k, Gprev);
      // S_work = rank_complete_s_old(S, Gprev, B.cols(0, k - 1));
      S_work = rank_complete_s(S, X, Gprev, B.cols(0, k - 1));
    }

    // prepare the Gamma slice to pass
    arma::mat Gprev_step;
    if (orth && k > 0) {
      Gprev_step = Gamma.cols(0, k - 1);
    } else {
      Gprev_step.reset();  // empty matrix
    }

    // one flip–flop update
    auto CAPre = CAP_one_component(S_work, X, T, beta_k, gamma_k, Gprev_step,
                                   max_iter, tol);

    // store results
    Gamma.col(k) = CAPre.gamma;
    B.col(k) = CAPre.beta;
    // copy back into Gamma_prev
    // Gamma_prev.col(k) = CAPre.gamma;
  }

  return Rcpp::List::create(Rcpp::Named("B") = B, Rcpp::Named("Gamma") = Gamma);
}

Rcpp::List CAP_multi_components(const arma::cube& S, const arma::mat& X,
                                const arma::vec& T, const int K,
                                const int max_iter = 200,
                                const double tol = 1e-6,
                                const bool orth = true) {
  const arma::uword p = S.n_rows;
  const arma::uword q = X.n_cols;

  arma::mat Gamma(p, K, arma::fill::zeros);  // store γ^{(1)},…,γ^{(K)}
  arma::mat B(q, K, arma::fill::zeros);      // store β^{(1)},…,β^{(K)}

  arma::mat Gamma_prev(p, K, arma::fill::zeros);  // pre-allocated

  arma::cube S_work = S;  // may be rank-completed each step

  for (int k = 0; k < K; ++k) {
    // initialize β and γ
    arma::vec beta_k(q, arma::fill::zeros);
    arma::vec gamma_k = arma::randn<arma::vec>(p);
    // gamma_k /= arma::norm(gamma_k, 2);

    // rank-complete S if orthogonality step required
    if (k > 0) {
      // use first k columns of Gamma_prev
      arma::mat Gprev = Gamma.cols(0, k - 1);
      gamma_k = orthogonalise_qr(gamma_k, Gprev);
      // S_work = rank_complete_s_old(S, Gprev, B.cols(0, k - 1));
      S_work = rank_complete_s(S, X, Gprev, B.cols(0, k - 1));
    }

    // prepare the Gamma slice to pass
    arma::mat Gprev_step;
    if (orth && k > 0) {
      Gprev_step = Gamma.cols(0, k - 1);
    } else {
      Gprev_step.reset();  // empty matrix
    }

    // one flip–flop update
    auto CAPre = CAP_one_component(S_work, X, T, beta_k, gamma_k, Gprev_step,
                                   max_iter, tol);

    // store results
    Gamma.col(k) = CAPre.gamma;
    B.col(k) = CAPre.beta;
    // copy back into Gamma_prev
    // Gamma_prev.col(k) = CAPre.gamma;
  }

  return Rcpp::List::create(Rcpp::Named("B") = B, Rcpp::Named("Gamma") = Gamma);
}
