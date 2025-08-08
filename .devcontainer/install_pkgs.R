## ----- non-interactive installs --------------------------------------------
## 1. Tell R where CRAN lives (Posit/RStudio’s “cloud” mirror is fast + HTTPS)
options(repos = c(CRAN = "https://cloud.r-project.org"))


target_lib <- .libPaths()[1]

## 3. Install in one shot, quietly, with dependencies resolved
install.packages(
    c(
        "Rcpp", "RcppArmadillo", "languageserver",
        "testthat", "roxygen2", "rversions", "urlchecker", 
        "styler", "xml2", "lintr", "devtools"
    ),
    quiet = FALSE, # suppress progress output
    dependencies = TRUE # pull in any needed supporting pkgs
)
