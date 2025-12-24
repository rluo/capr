library(devtools)

devtools::document()
devtools::install()

devtools::test()

library(capr)


simu.data <- simu.capr(seed = 123L, n = 120L)
K <- 4L
fit <- capr(
    S = simu.data$S,
    X = simu.data$X,
    K = K,
    orth = TRUE
)
fit



plot(fit, simu.data$S)
