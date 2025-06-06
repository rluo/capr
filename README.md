# capa

R package for implementing the algorithms in the CAP [paper](https://doi.org/10.1093/biostatistics/kxz057) using [Armadillo](https://arma.sourceforge.net/).



# Development environment

Package skeleton:

    library(usethis)
    create_package(".")

Add Rcpp support:

    usethis::use_rcpp("capa")
    usethis::use_package("RcppArmadillo", type = "LinkingTo")
