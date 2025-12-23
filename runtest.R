library(devtools)

devtools::document()
devtools::install()

devtools::test()

library(capr)
simu.data <- simu.capr(seed = 123L, p = 5L, n = 120L)
K <- 2L

print("No orth testing ")


fit <- capr(
    S = simu.data$S,
    X = simu.data$X,
    K = K,
    weight = runif(simu.data$n) + 0.5,
    orth = FALSE
)
fit


print("With orth testing ")

fit <- capr(
    S = simu.data$S,
    X = simu.data$X,
    K = K,
    weight = runif(simu.data$n) + 0.5,
    orth = TRUE
)
fit
