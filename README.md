# capa

R package for implementing the algorithms in the CAP [paper](https://doi.org/10.1093/biostatistics/kxz057) using [Armadillo](https://arma.sourceforge.net/).



# Development environment

Package skeleton:

    library(usethis)
    create_package(".")
    usethis::create_package("path/to/myrcpppkg")

Add Rcpp support:

    usethis::use_rcpp("capr")
    usethis::use_rcpp()                                       # sets up Rcpp glue
    usethis::use_package("Rcpp", type = "LinkingTo")
    usethis::use_package("Rcpp", type = "Imports")
    usethis::use_package("RcppArmadillo", type = "LinkingTo")

Add Roxygen

    usethis::use_roxygen_md()


Add Testing

    usethis::use_testthat(edition = 3)

Add lintr

    usethis::use_lintr()
