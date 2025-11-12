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

/**
 * @title Rank-completion of Sigma matrices via the max-exp rule
 * @name rank_complete_s
 * @description
 * Rank-completes a stack of covariance (Sigma) matrices using the
 * “max-exp” rule from a generalized low-rank structure. Given
 * previously-fitted directions \eqn{\Gamma} and coefficients \eqn{B},
 * it computes \eqn{E = \exp(X B)} (elementwise), takes
 * \eqn{a = \max_{i,k} E_{ik}}, and for each slice \eqn{i=1,\dots,n}
 * returns
 * \deqn{\tilde{S}_i \;=\; S_i \;+\; \Gamma \; \mathrm{diag}(a - E_{i\cdot}) \;
 * \Gamma^\top.}
 *
 * @details
 * **Inputs and shapes (column-major, Armadillo default)**
 * - \code{S_cube}: \eqn{p \times p \times n} stack of original
 *   covariance matrices \eqn{S_i} (symmetric positive definite).
 * - \code{X}: \eqn{n \times q} design matrix (rows are \eqn{x_i^\top}).
 * - \code{Gamma_prev}: \eqn{p \times (k-1)} matrix of previously fitted
 *   directions \eqn{\Gamma^{(k-1)}}.
 * - \code{B}: \eqn{q \times (k-1)} coefficient matrix whose
 *   \eqn{k-1} columns correspond to components; its first row may
 *   represent an intercept if your model uses one (then ensure
 *   \code{X} contains a column of ones).
 *
 * **Algorithm (max-exp rule)**
 * 1. \eqn{E \leftarrow \exp(X B)} (elementwise), \eqn{E \in \mathbb{R}^{n
 * \times (k-1)}}.
 * 2. \eqn{a \leftarrow \max(E)} (scalar maximum over all entries).
 * 3. For each \eqn{i=1,\dots,n}, set \eqn{d_i = a \mathbf{1} - E_{i\cdot}} and
 *    \eqn{\tilde{S}_i = S_i + \Gamma \mathrm{diag}(d_i)\Gamma^\top}.
 *
 * **Edge cases**
 * - If \eqn{k-1 = 0} (i.e., \code{Gamma_prev} has zero columns),
 *   the input \code{S_cube} is returned unchanged.
 *
 * **Dimension checks**
 * - \code{nrow(X)} must equal the number of slices \code{n} in \code{S_cube}.
 * - \code{B} must be \eqn{q \times (k-1)} where \eqn{q = ncol(X)}.
 *
 * @param S_cube A numeric cube of size \eqn{p \times p \times n}:
 *   original covariance matrices \eqn{S_i}.
 * @param X A numeric matrix \eqn{n \times q}: design matrix (rows
 * \eqn{x_i^\top}).
 * @param Gamma_prev A numeric matrix \eqn{p \times (k-1)}:
 *   previously fitted directions \eqn{\Gamma^{(k-1)}}.
 * @param B A numeric matrix \eqn{q \times (k-1)}: coefficient matrix
 *   (column \eqn{j} corresponds to component \eqn{j}).
 *
 * @return A numeric cube \eqn{p \times p \times n} containing the completed
 * matrices \eqn{\tilde{S}_i^{(k)}} in the same slice order as \code{S_cube}.
 *
 * @section Notes:
 * - All matrices are treated as column-major.
 * - No symmetry-enforcement step is applied; inputs should be symmetric
 *   and operations preserve symmetry up to numerical rounding.
 *
 * @examples
 * \dontrun{
 *   set.seed(1)
 *   p <- 3; n <- 4; q <- 2; k1 <- 2  # k-1 = 2
 *   S_cube <- array(0, dim = c(p, p, n))
 *   for (i in 1:n) {
 *     A <- matrix(rnorm(p*p), p, p)
 *     S_cube[,,i] <- crossprod(A) + diag(p)  # SPD
 *   }
 *   X <- matrix(rnorm(n*q), n, q)
 *   Gamma_prev <- matrix(rnorm(p*k1), p, k1)
 *   B <- matrix(rnorm(q*k1), q, k1)
 *   out <- rank_complete_s(S_cube, X, Gamma_prev, B)
 *   dim(out)  # p x p x n
 * }
 *
 * @seealso \code{\link[=Rcpp]{Rcpp}}, \code{\link{RcppArmadillo}}
 * @export
 */
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

/**
 * @title Cosine Similarity Between Two Vectors
 * @name cosine_similarity
 * @description
 * Computes the cosine similarity between two numeric vectors
 * \eqn{a} and \eqn{b}:
 * \deqn{\cos(\theta) \;=\; \frac{a^\top b}{\|a\|_2 \; \|b\|_2}}
 * where \eqn{\|\cdot\|_2} denotes the Euclidean (L2) norm.
 *
 * @details
 * - The function accepts any Armadillo dense vector type (e.g. `arma::vec`,
 *   `arma::rowvec`, or subviews).
 * - Both vectors must be the same length; otherwise an error is thrown.
 * - If either vector has (near-)zero norm (smaller than \code{eps}),
 *   the function throws an error.
 *
 * @param a A numeric vector.
 * @param b A numeric vector of the same length as \code{a}.
 * @param eps A small numeric tolerance (default \code{1e-12}) used to guard
 *   against division by zero if one of the vectors is near zero length.
 *
 * @return A scalar numeric value in \eqn{[-1, 1]} giving the cosine similarity.
 *
 * @examples
 * # Two identical vectors => cosine similarity = 1
 * cosine_similarity(c(1, 2, 3), c(1, 2, 3))
 *
 * # Orthogonal vectors => cosine similarity = 0
 * cosine_similarity(c(1, 0), c(0, 1))
 *
 * # Opposite vectors => cosine similarity = -1
 * cosine_similarity(c(1, 2), c(-1, -2))
 *
 * @seealso \code{\link[=Rcpp]{Rcpp}}, \code{\link{RcppArmadillo}}
 * @export
 */
// [[Rcpp::export]]
double cosine_similarity(const arma::vec& a, const arma::vec& b,
                         const double eps = 1e-12) {
  if (a.n_elem != b.n_elem) Rcpp::stop("cosine_similarity: size mismatch");

  const double na = arma::norm(a, 2);
  const double nb = arma::norm(b, 2);

  if (na < eps || nb < eps) Rcpp::stop("cosine_similarity: zero vector");

  return arma::dot(a, b) / (na * nb);
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
