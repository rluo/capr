# capr

R package for implementing the algorithms in the CAP [paper](https://doi.org/10.1093/biostatistics/kxz057) using [Armadillo](https://arma.sourceforge.net/).


# Project Roadmap

## Current Focus (Version 0.2)

*  Implement core functionality using Armadillo and C++ 17
   *  [ ] Core CAP regression with and without orthogonality constraints
   *  [ ] S3 print and summary mehtods
   *  [ ] Improved component effects removal method


## Future Releases

*  High dimensional CAP
*  Longitudinal CAP



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
