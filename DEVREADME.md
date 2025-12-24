

# Package build

    devtools::document()
    devtools::test()
    devtools::install()
    devtools::check()


# Development environment setup

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


Add Releases and TODOs

    usethis::use_github()
    usethis::use_tidy_issue_template()
    usethis::use_news_md()


# Test

    Rscript --vanilla  runtest.R

## Release

    usethis::use_version("minor")   # or patch/major
    usethis::use_cran_comments()
    usethis::use_release_issue()    # checklist on GitHub
    devtools::check()
    # After tagging and pushing:
    usethis::use_github_release()