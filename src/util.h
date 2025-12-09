#pragma once
#include <RcppArmadillo.h>

double cosine_similarity_cpp(const arma::vec& a, const arma::vec& b,
                             double eps = 1e-12);

double log_deviation_from_diagonality_cpp(const arma::cube& S_cube,
                                          const arma::vec& nval,
                                          const arma::mat& B);

double cap_loglike_cpp(const arma::cube& S_cube, const arma::mat& X,
                       const arma::vec& T, const arma::vec& beta,
                       const arma::vec& gamma);