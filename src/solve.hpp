#pragma once
#include <RcppArmadillo.h>

//--------------------------------------------------------------------------
//  γ̂  for the unconstrained case  (Γ_prev has zero columns)
//
//    minimise   γᵀ A γ      s.t.   γᵀ H γ = 1
//
//  •  H  must be symmetric positive-definite
//  •  A  symmetric
//--------------------------------------------------------------------------
static arma::vec solve_gamma_unconstrained(const arma::mat& A,
                                           const arma::mat& H) {
  // --- spectral square-root of H  ----------------------------------------
  arma::vec hval;
  arma::mat Hevec;
  arma::eig_sym(hval, Hevec, H, "std");  // H = V Λ Vᵀ  (Λ>0)

  arma::mat H_half = Hevec * arma::diagmat(arma::sqrt(hval)) * Hevec.t();
  arma::mat H_invhalf =
      Hevec * arma::diagmat(1.0 / arma::sqrt(hval)) * Hevec.t();

  // --- ordinary eigenproblem  B ν = λ ν  ---------------------------------
  arma::mat B = H_invhalf * A * H_invhalf;  // B symmetric
  arma::vec lam;
  arma::mat V;
  arma::eig_sym(lam, V, B);     // ascending eigenvalues
  arma::vec nu_min = V.col(0);  // ν with smallest λ

  // --- transform back   γ = H⁻½ ν   and normalise  γᵀHγ = 1  --------------
  arma::vec gamma = H_invhalf * nu_min;
  gamma /= std::sqrt(arma::as_scalar(gamma.t() * H * gamma));

  return gamma;
}

// ─────────────────────────────────────────────────────────────────────────────
//  γ–step  :  minimise  γᵀAγ   s.t.  γᵀHγ = 1  and  Γᵀγ = 0               (*)
//  when Γ is empty: fall back to the unconstrained solver
// ─────────────────────────────────────────────────────────────────────────────
static arma::vec solve_gamma(
    const arma::mat& A,           //  p × p   (symmetric, SPD-ish)
    const arma::mat& H,           //  p × p   (SPD)
    const arma::mat& Gamma_prev)  //  p × (k-1)  (may have 0 cols)
{
  const arma::uword p = A.n_rows;

  // -------------------------------------------------------------------------
  //  Case 1  :  *no* orthogonality constraints  →  plain gen-eigen problem
  // -------------------------------------------------------------------------
  //   if (Gamma_prev.n_cols == 0) {
  //     // ← new unconstrained branch
  //     return solve_gamma_unconstrained(A, H);
  //   }

  // -------------------------------------------------------------------------
  //  Case 2  :  orthogonality w.r.t. Γ_prev  (same code as before)
  // -------------------------------------------------------------------------
  arma::vec eval;
  arma::mat evec;
  arma::eig_sym(eval, evec, H, "std");  // H = V Λ Vᵀ

  arma::mat H_half = evec * arma::diagmat(arma::sqrt(eval)) * evec.t();
  arma::mat H_invhalf = evec * arma::diagmat(1.0 / arma::sqrt(eval)) * evec.t();

  arma::mat B = H_invhalf * A * H_invhalf;  // H⁻½ A H⁻½
  arma::mat U = H_invhalf * Gamma_prev;     // p × (k-1)

  arma::mat Q = arma::null(U.t());  // p × (p-k+1)

  arma::mat Bred = Q.t() * B * Q;  // reduced matrix
  arma::vec lam;
  arma::mat Z;
  arma::eig_sym(lam, Z, Bred, "std");
  arma::vec z_min = Z.col(0);  // smallest eigen-pair

  arma::vec gamma = H_invhalf * Q * z_min;
  gamma /= std::sqrt(arma::as_scalar(gamma.t() * H * gamma));  // γᵀHγ = 1
  return gamma;
}

// ─────────────────────────────────────────────────────────────────────────────
//  β–STEP  (one Newton / Fisher like scoring update)
// ─────────────────────────────────────────────────────────────────────────────
static arma::vec newton_beta(const arma::cube& S, const arma::mat& X,
                             const arma::vec& T, const arma::vec& beta,
                             const arma::vec& gamma) {
  const arma::uword n = S.n_slices;
  arma::mat Hbeta(beta.n_elem, beta.n_elem, arma::fill::zeros);
  arma::vec g(beta.n_elem, arma::fill::zeros);

  for (arma::uword i = 0; i < n; ++i) {
    double eta = arma::dot(X.row(i), beta);
    double wi =
        std::exp(-eta) * arma::as_scalar(gamma.t() * S.slice(i) * gamma);
    arma::vec xi = X.row(i).t();
    Hbeta += wi * (xi * xi.t());
    g += (T[i] - wi) * xi;
  }
  return beta - arma::solve(Hbeta, g, arma::solve_opts::likely_sympd);
}
