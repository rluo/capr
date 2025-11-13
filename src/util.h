#pragma once
#include <RcppArmadillo.h>

double cosine_similarity(const arma::vec& a, const arma::vec& b,
                         double eps = 1e-12);
