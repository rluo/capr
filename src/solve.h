#pragma once
#include <RcppArmadillo.h>

inline arma::mat make_symmetric(const arma::mat& M) {
  return 0.5 * (M + M.t());
}

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
  arma::mat Hsym = make_symmetric(H);
  arma::vec hval;
  arma::mat Hevec;
  arma::eig_sym(hval, Hevec, Hsym, "std");  // H = V Λ Vᵀ  (Λ>0)

  arma::mat H_half = Hevec * arma::diagmat(arma::sqrt(hval)) * Hevec.t();
  arma::mat H_invhalf =
      Hevec * arma::diagmat(1.0 / arma::sqrt(hval)) * Hevec.t();

  // --- ordinary eigenproblem  B ν = λ ν  ---------------------------------
  arma::mat Asym = make_symmetric(A);
  arma::mat B = make_symmetric(H_invhalf * Asym * H_invhalf);  // B symmetric
  arma::vec lam;
  arma::mat V;
  arma::eig_sym(lam, V, B);     // ascending eigenvalues
  arma::vec nu_min = V.col(0);  // ν with smallest λ

  // --- transform back   γ = H⁻½ ν   and normalise  γᵀHγ = 1  --------------
  arma::vec gamma = H_invhalf * nu_min;
  gamma /= std::sqrt(arma::as_scalar(gamma.t() * Hsym * gamma));

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
  // const arma::uword p = A.n_rows;

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
  arma::mat Hsym = make_symmetric(H);
  arma::vec eval;
  arma::mat evec;
  arma::eig_sym(eval, evec, Hsym, "std");  // H = V Λ Vᵀ

  arma::mat H_half = evec * arma::diagmat(arma::sqrt(eval)) * evec.t();
  arma::mat H_invhalf = evec * arma::diagmat(1.0 / arma::sqrt(eval)) * evec.t();

  arma::mat Asym = make_symmetric(A);
  arma::mat B = make_symmetric(H_invhalf * Asym * H_invhalf);  // H⁻½ A H⁻½
  arma::mat U = H_invhalf * Gamma_prev;                        // p × (k-1)

  arma::mat Q = arma::null(U.t());  // p × (p-k+1)

  arma::mat Bred = make_symmetric(Q.t() * B * Q);  // reduced matrix
  arma::vec lam;
  arma::mat Z;
  arma::eig_sym(lam, Z, Bred, "std");
  arma::vec z_min = Z.col(0);  // smallest eigen-pair

  arma::vec gamma = H_invhalf * Q * z_min;
  gamma /= std::sqrt(arma::as_scalar(gamma.t() * Hsym * gamma));  // γᵀHγ = 1
  return gamma;
}

// ─────────────────────────────────────────────────────────────────────────────
//  β–STEP  (one Newton / Fisher like scoring update)
// ─────────────────────────────────────────────────────────────────────────────
arma::vec newton_beta(const arma::cube& S, const arma::mat& X,
                      const arma::vec& T, const arma::vec& beta_init,
                      const arma::vec& gamma, int max_iter = 1,
                      double tol = 1e-6);
