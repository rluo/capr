#pragma once
#include <RcppArmadillo.h>

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
  const double a = 2 * E.max();    // scalar

  // ---- pre-allocate result cube ------------------------------------------
  arma::cube S_tilde(p, p, n, arma::fill::zeros);

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

arma::cube deflate_s_projection(const arma::cube& S_cube,
                                const arma::mat& Gamma_prev) {
  const arma::uword p = S_cube.n_rows;
  const arma::uword n = S_cube.n_slices;

  if (Gamma_prev.n_cols == 0) return S_cube;

  arma::mat Q, R;
  arma::qr_econ(Q, R, Gamma_prev);  // Q: p×(k−1), orthonormal columns
  arma::mat P = Q * Q.t();          // true orthogonal projector
  arma::mat I = arma::eye(p, p);
  arma::mat M = I - P;

  arma::cube S_tilde(p, p, n, arma::fill::zeros);
  for (arma::uword i = 0; i < n; ++i) {
    arma::mat Si = S_cube.slice(i);
    S_tilde.slice(i) = M * Si * M;  // (I-P) S (I-P)
    S_tilde.slice(i) = 0.5 * (S_tilde.slice(i) + S_tilde.slice(i).t());  // sym
  }
  return S_tilde;
}

arma::mat inv_sqrt_matrix(const arma::mat& A, const int k) {
  // Force symmetry before the eigendecomposition to avoid warnings/failures
  arma::mat A_sym = 0.5 * (A + A.t());

  arma::vec eigval;
  arma::mat eigvec;
  const bool ok = arma::eig_sym(eigval, eigvec, A_sym, "dc");
  if (!ok) Rcpp::stop("inv_sqrt_matrix: eigen-decomposition failed.");

  // Guard against non-positive eigenvalues (rank deficiency / numerical noise)
  const double lambda_max = eigval.max();
  const double floor =
      std::max(lambda_max * 1e-6, 1e-12);  // cap condition number, avoid zeros
  eigval.transform([floor](double v) { return (v < floor) ? floor : v; });

  arma::vec inv_sqrt_eigval = 1.0 / arma::sqrt(eigval.head(k));
  arma::mat inv_sqrt_D = arma::diagmat(inv_sqrt_eigval);
  arma::mat eigvec_k = eigvec.cols(0, k - 1);

  return eigvec_k * inv_sqrt_D * eigvec_k.t();
}

arma::cube rank_complete_multiply(const arma::cube& S_cube,     // p×p×n
                                  const arma::mat& X,           // n×q
                                  const arma::mat& Gamma_prev,  // p×(k−1)
                                  const arma::mat& B)           // q×(k−1)
{
  const arma::uword n = S_cube.n_slices;
  const arma::uword p = S_cube.n_rows;

  if (Gamma_prev.n_cols == 0) return S_cube;

  if (X.n_rows != n) Rcpp::stop("rank_complete_multiply: X must have n rows.");
  if (B.n_rows != X.n_cols || B.n_cols != Gamma_prev.n_cols)
    Rcpp::stop(
        "rank_complete_multiply: B must be q×(k-1) and match Gamma_prev.");

  arma::mat E = arma::exp(X * B);  // n × (k−1)

  // ---- pre-allocate result cube ------------------------------------------
  arma::cube S_tilde(p, p, n, arma::fill::zeros);

  // ---- main loop over subjects -------------------------------------------
  for (arma::uword i = 0; i < n; ++i) {
    arma::rowvec Ei = E.row(i);  // 1 × (k−1)

    arma::mat Shat = Gamma_prev * arma::diagmat(Ei) * Gamma_prev.t();  // p × p
    // Shat = 0.5 * (Shat + Shat.t());  // enforce symmetry for eig_sym
    arma::mat Ssqrtinv = inv_sqrt_matrix(Shat, Gamma_prev.n_cols);

    S_tilde.slice(i) = Ssqrtinv * S_cube.slice(i) * Ssqrtinv;
  }

  return S_tilde;
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
