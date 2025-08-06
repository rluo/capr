library(RcppArmadillo)

sourceCpp("dev/cap_c.cpp")
source("cap_c.R")
## (Assuming the two files are in your package or have been sourced)
n <- 20
p <- 3
Ti <- 50
set.seed(123)
X <- cbind(1, rnorm(n)) # intercept + one covariate
Y <- lapply(seq_len(n), function(i) matrix(rnorm(Ti * p), Ti, p))

out <- cap_first_direction(Y, X)
str(out)
