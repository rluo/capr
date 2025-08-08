// [[Rcpp::depends(RcppArmadillo)]]
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
  if (Gamma_prev.n_cols == 0) {
    // ← new unconstrained branch
    return solve_gamma_unconstrained(A, H);
  }

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
//  β–STEP  (one Newton / Fisher scoring update)
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

// ─────────────────────────────────────────────────────────────────────────────
//  Thin POD used for return-by-value  (RVO friendly)
// ─────────────────────────────────────────────────────────────────────────────
struct CAPResult {
  arma::vec beta;
  arma::vec gamma;
};

// ─────────────────────────────────────────────────────────────────────────────
//   One CAP component – flip-flop,  now returning  CAPResult
// ─────────────────────────────────────────────────────────────────────────────
static CAPResult CAP_one_component(const arma::cube& S, const arma::mat& X,
                                   const arma::vec& T,
                                   const arma::vec& beta_init,
                                   const arma::vec& gamma_init,
                                   const arma::mat& Gamma_prev,
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

  for (int it = 0; it < max_iter; ++it) {
    arma::vec beta_old = beta, gamma_old = gamma;

    beta = newton_beta(S, X, T, beta, gamma);

    // build A(β)
    arma::mat A(p, p, arma::fill::zeros);
    for (arma::uword i = 0; i < n; ++i)
      A += std::exp(-arma::dot(X.row(i), beta)) * S.slice(i);

    gamma = solve_gamma(A, H, Gamma_prev);

    if (std::max(arma::norm(beta - beta_old, "inf"),
                 arma::norm(gamma - gamma_old, "inf")) < tol)
      break;
  }

  return {std::move(beta), std::move(gamma)};  // NRVO
}

/*──────────────────────────────────────────────────────────────────────────────
  Rank–completion step used in the CAP flip-flop:

      Ŷᵢ   …  p×p data matrices   (cube  p×p×n)
      Γ     …  p×(k-1)  previously-estimated directions
      B     …  q×(k-1)  coefficient matrix,   row 0 = intercepts β_{j0}

  For each i = 1..n
  ────────────────────────────────────────────────────────────────────────────
    1.   Ŷᵢ  = Yᵢ  −  Yᵢ Γ Γᵀ                       (residual)
    2.   U D Vᵀ = svd( Ŷᵢ )
    3.   D̃  = diag(  D₁ … D_{p-(k-1)},  √e^{β_{10}}, …, √e^{β_{(k-1)0}} )
    4.   Ỹᵢ = U D̃ Vᵀ
  ────────────────────────────────────────────────────────────────────────────
  Returns a cube of the same shape as Y with slices Ỹᵢ.
────────────────────────────────────────────────────────────────────────────*/

// [[Rcpp::export]]
arma::cube rank_complete_y(const arma::cube& Y,          // p × p × n
                           const arma::mat& Gamma_prev,  // p × (k-1)
                           const arma::mat& B)           // q × (k-1)
{
  const arma::uword p = Y.n_rows;
  const arma::uword n = Y.n_slices;
  const arma::uword k_minus = Gamma_prev.n_cols;  // = k−1

  if (k_minus != 0 && B.n_cols != k_minus)
    Rcpp::stop("B must have k-1 columns, matching Gamma_prev.");

  arma::cube Y_tilde(p, p, n);

  // intercepts β_{j0}  (first row of B)
  arma::rowvec beta0 = (k_minus == 0) ? arma::rowvec() : B.row(0);

  for (arma::uword i = 0; i < n; ++i) {
    arma::mat Yi = Y.slice(i);

    // 1. residual   Ŷᵢ
    arma::mat Yhat =
        (k_minus == 0) ? Yi : Yi - Yi * Gamma_prev * Gamma_prev.t();

    // 2. thin SVD   Ŷᵢ = U D Vᵀ
    arma::mat U, V;
    arma::vec D;
    arma::svd(U, D, V, Yhat, "std");

    // 3. build D̃
    arma::vec Dtilde(p, arma::fill::zeros);
    arma::uword keep = p - k_minus;  // #original sing.vals

    if (keep > 0) Dtilde.head(keep) = D.head(keep);

    for (arma::uword j = 0; j < k_minus; ++j)
      Dtilde(keep + j) = std::sqrt(std::exp(beta0(j)));

    // 4. reconstruct
    Y_tilde.slice(i) = U * arma::diagmat(Dtilde) * V.t();
  }

  return Y_tilde;
}

/*──────────────────────────────────────────────────────────────────────────────
   ❶  Rank–completion for the Σ-matrices  (Section 5.2, “fix S_i”)

     •  Yields  S̃_i^(k)  that carry over variance in the Γ^(k-1) directions
        via   exp(β_{j0})  (intercepts in B).

   INPUT
     S_cube      p × p × n    original S_i  (symmetric, p.d.)
     Gamma_prev  p × (k-1)    previously fitted directions
     B           q × (k-1)    coefficient matrix, row-0 = β_{j0}

   OUTPUT
     p × p × n   cube S_tilde  (same shape as S_cube)

   ALGORITHM  (for each i = 1..n)
     P      = Γ Γᵀ
     S_res  = (I−P) S_i (I−P)                   # rank ≤ p−(k−1)
     [V,D]  = eig_sym(S_res)                    # D ascending
     Replace the first (k−1) eigen-values       # (smallest / ≈ 0)
         D̃_j = exp( β_{j0} )                   #  j = 1..k−1
     S̃_i   = V diag(D̃) Vᵀ
──────────────────────────────────────────────────────────────────────────────*/

// [[Rcpp::export]]
arma::cube rank_complete_s_old(const arma::cube& S_cube,     // p×p×n
                               const arma::mat& Gamma_prev,  // p×(k−1)
                               const arma::mat& B)           // q×(k−1)
{
  const arma::uword p = S_cube.n_rows, n = S_cube.n_slices,
                    kminus1 = Gamma_prev.n_cols;  // k-1

  if (kminus1 != 0 && B.n_cols != kminus1)
    Rcpp::stop("B must have k-1 columns, same as Gamma_prev.");

  // intercepts β_{j0}  (first row)
  arma::rowvec beta0 = (kminus1 == 0) ? arma::rowvec() : B.row(0);

  const arma::mat P = (kminus1 == 0) ? arma::mat(p, p, arma::fill::zeros)
                                     : Gamma_prev * Gamma_prev.t();
  const arma::mat IminusP = arma::eye(p, p) - P;

  arma::cube S_tilde(p, p, n);

  for (arma::uword i = 0; i < n; ++i) {
    arma::mat S_res = IminusP * S_cube.slice(i) * IminusP;  // (1)

    arma::vec eval;
    arma::mat evec;
    arma::eig_sym(eval, evec, S_res, "std");  // ascending λ₁ ≤ … ≤ λ_p

    // Replace the k−1 smallest eigenvalues
    for (arma::uword j = 0; j < kminus1; ++j)
      eval(j) = std::exp(beta0(j));  //   √ dropped in derivation

    arma::mat S_new = evec * arma::diagmat(eval) * evec.t();
    S_tilde.slice(i) =
        0.5 * (S_new + S_new.t());  // enforce symmetry numerically
  }
  return S_tilde;
}

//////////////////////////////////////////////////////////////////////////////
//  Rank-completion of Σ-matrices via “max-exp” rule (latex description)
//
//  INPUT
//    • S_cube      p × p × n      original   Σ_i      (symmetric p.d.)
//    • X           n × q          design matrix      (rows x_i^T)
//    • Gamma_prev  p × (k−1)      previously fitted directions Γ^{(k-1)}
//    • B           q × (k−1)      coefficient matrix (col-j for component j)
//                                  — first row is β_{j0} (intercept)
//
//  ALGORITHM
//      E_{ik}  = exp( x_i^T  B_{·k} )                         // n × (k−1)
//      a       = max_{i,k}  E_{ik}
//
//      for i = 1..n
//          d_i  = a − E_{i·}          // length (k−1)
//          S̃_i = S_i  +  Γ diag(d_i) Γᵀ
//
//  OUTPUT
//      cube S_tilde   p × p × n       completed Σ̃_i^(k)
//
//  NOTES
//      • If k−1 = 0 (i.e. Γ_prev has zero columns) the function simply
//        returns S_cube unchanged.
//      • All matrices are treated as column-major (Armadillo default).
//////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::cube rank_complete_s(const arma::cube& S_cube,     // p×p×n
                           const arma::mat& X,           // n×q
                           const arma::mat& Gamma_prev,  // p×(k−1)
                           const arma::mat& B)           // q×(k−1)
{
  const arma::uword p = S_cube.n_rows;
  const arma::uword n = S_cube.n_slices;
  const arma::uword k_minus1 = Gamma_prev.n_cols;  // k−1

  // ---- trivial case : no previous directions -----------------------------
  if (k_minus1 == 0) return S_cube;

  // ---- dimension checks ---------------------------------------------------
  if (X.n_rows != n) Rcpp::stop("rank_complete_s_maxexp: X must have n rows.");
  if (B.n_rows != X.n_cols || B.n_cols != k_minus1)
    Rcpp::stop("rank_complete_s_maxexp: B must be q×(k-1) and match X.");

  // ---- compute   E = exp(X %*% B)   and  a = max(E) -----------------------
  arma::mat E = arma::exp(X * B);  // n × (k−1)
  const double a = E.max();        // scalar

  // ---- pre-allocate result cube ------------------------------------------
  arma::cube S_tilde(p, p, n);

  // ---- main loop over subjects -------------------------------------------
  for (arma::uword i = 0; i < n; ++i) {
    arma::rowvec Ei = E.row(i);   // 1 × (k−1)
    arma::vec di = (a - Ei).t();  // (k−1) × 1

    arma::mat correction =
        Gamma_prev * arma::diagmat(di) * Gamma_prev.t();  // p × p

    S_tilde.slice(i) = S_cube.slice(i) + correction;
  }

  return S_tilde;
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
Rcpp::List cap_one_comp(const arma::cube& S, const arma::mat& X,
                        const arma::vec& T, arma::vec beta0, arma::vec gamma0,
                        const arma::mat& Gamma_prev, int max_iter = 200,
                        double tol = 1e-6) {
  CAPResult res =
      CAP_one_component(S, X, T, beta0, gamma0, Gamma_prev, max_iter, tol);

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

//--------------------------------------------------------------------------
//  Cosine similarity   cos(θ) = (aᵀ b) / (‖a‖₂ ‖b‖₂)
//  • Accepts any Armadillo dense vector flavour (vec / rowvec / subview)
//  • Throws if one of the inputs has (near-)zero length.
//--------------------------------------------------------------------------
// [[Rcpp::export]]
double cosine_similarity(const arma::vec& a, const arma::vec& b,
                         const double eps = 1e-12) {
  if (a.n_elem != b.n_elem) Rcpp::stop("cosine_similarity: size mismatch");

  const double na = arma::norm(a, 2);
  const double nb = arma::norm(b, 2);

  if (na < eps || nb < eps) Rcpp::stop("cosine_similarity: zero vector");

  return arma::dot(a, b) / (na * nb);
}
// [[Rcpp::export]]
Rcpp::List cap_multi_comp(const arma::cube& S, const arma::mat& X,
                          const arma::vec& T, const int K,
                          const int max_iter = 200, const double tol = 1e-6,
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
