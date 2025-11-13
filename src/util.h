#pragma once
#include <RcppArmadillo.h>

double cosine_similarity(const arma::vec& a, const arma::vec& b,
                         double eps = 1e-12);

double log_deviation_from_diagonality(const arma::cube& S_cube,
                                      const arma::vec& nval,
                                      const arma::mat& B);
