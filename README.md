# capp

R package for implementing the algorithms in the CAP [paper](https://doi.org/10.1093/biostatistics/kxz057)



# Development environment

Package skeleton:

    library(usethis)
    create_package(".")

Add Rcpp support:

    usethis::use_rcpp("capp")
    usethis::use_package("RcppArmadillo", type = "LinkingTo")
