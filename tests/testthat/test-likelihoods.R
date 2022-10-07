library(testthat)
library(EmpiricalCalibration)


# test Gaussian product
test_that("gaussianProduct - input not numeric", {
  expect_error(gaussianProduct("a", "b", "c", "d"))
  expect_equal(gaussianProduct(1, 2, 0.1, NaN), NaN)
  expect_equal(gaussianProduct(1, 2, 0.1, Inf), 0)
})

test_that("gaussianProduct - both std are zero", {
  mu1 <- 1
  mu2 <- 2
  sd1 <- 0
  sd2 <- 0
  result <- gaussianProduct(mu1, mu2, sd1, sd2)
  expect_equal(result, NaN)
})
test_that("gaussianProduct result is correct", {
  mu1 <- 1
  mu2 <- 2
  sd1 <- 0.1
  sd2 <- 0.01
  result <- gaussianProduct(mu1, mu2, sd1, sd2)
  expect_equal(result, 1.256095e-21)
})

# test logLikelihoodNull
test_that("logLikelihoodNull - negative theta2", {
  result <- EmpiricalCalibration:::logLikelihoodNull(theta = c(1, -1), logRr = c(1, 2), seLogRr = c(10.2))
  expect_equal(result, Inf)
})
test_that("logLikelihoodNull - logRr and seLogRr input", {
  expect_equal(EmpiricalCalibration:::logLikelihoodNull(theta = c(1, 2), logRr = c(NaN, 2), seLogRr = c(1, 0.1)), NaN)
  expect_equal(EmpiricalCalibration:::logLikelihoodNull(theta = c(1, 2), logRr = c(1, 2), seLogRr = c(1, NaN)), NaN)
  expect_equal(EmpiricalCalibration:::logLikelihoodNull(theta = c(1, 2), logRr = c(Inf, 2), seLogRr = c(1, 0.1)), Inf)
  expect_equal(EmpiricalCalibration:::logLikelihoodNull(theta = c(1, 2), logRr = c(1, 2), seLogRr = c(1, Inf)), Inf)
})
test_that("logLikelihoodNull - correct calculation", {
  sd <- 1e-7
  result <- logLikelihoodNull(theta = c(1, 1 / sd^2), logRr = c(1, 2), seLogRr = c(0.01, 0.1))
  expect_equal(result, (0 - dnorm(1, 1, 0.01, log = TRUE) - dnorm(1, 2, 0.1, log = TRUE)))
})
test_that("logLikelihoodNull - result is infinite", {
  result <- logLikelihoodNull(theta = c(1, 0.001), logRr = c(1, 0), seLogRr = c(0.1, 2^1023))
  expect_equal(result, Inf)
})

# test logLikelihoodNullMcmc
test_that("logLikelihoodNullMcmc - logRr and seLogRr input", {
  expect_equal(logLikelihoodNullMcmc(theta = c(1, 2), logRr = c(NaN, 2), seLogRr = c(1, 0.1)), NaN)
  expect_equal(logLikelihoodNullMcmc(theta = c(1, 2), logRr = c(1, 2), seLogRr = c(1, NaN)), NaN)
})
# calculations are straightforward
test_that("logLikelihoodNullMcmc - result", {
  theta <- c(1, 2)
  logRr <- c(1, 2)
  seLogRr <- c(0.01, 0.1)
  result1 <- (0 - log(gaussianProduct(1, 1, 0.01, 1 / sqrt(2))) - log(gaussianProduct(2, 1, 0.1, 1 / sqrt(2))) - dgamma(2, 1e-04, 1e-04, log = TRUE))
  result <- logLikelihoodNullMcmc(theta, logRr, seLogRr)
  expect_equal(result, result1)
})

# test minLogLikelihoodErrorModel
test_that("minLogLikelihoodErrorModel - estimateLl function cases", {
  theta <- c(3, 4, -2, 5)
  logRr <- c(2)
  seLogRr <- c(0.2)
  trueLogRr <- c(0.25)
  expect_equal(EmpiricalCalibration:::minLogLikelihoodErrorModel(theta, logRr, seLogRr, trueLogRr), 99999)
  theta <- c(3, 4, 2e-7, 5e-7)
  logRr <- c(2)
  seLogRr <- c(0.2)
  trueLogRr <- c(0.25)
  expect_equal(minLogLikelihoodErrorModel(theta, logRr, seLogRr, trueLogRr), -dnorm(2, 4, 0.2, log = TRUE))
  theta <- c(3, 4, 2e-5, 5e-5)
  logRr <- c(2)
  seLogRr <- c(0.2)
  trueLogRr <- c(0.25)
  expect_equal(
    minLogLikelihoodErrorModel(theta, logRr, seLogRr, trueLogRr),
    -log(gaussianProduct(2, 4, 0.2, 3.25e-05))
  )
})

# minLogLikelihoodErrorModelLegacy
test_that("minLogLikelihoodErrorModelLegacy - is infinite result", {
  theta <- c(3, 4, 2e-5, 5e-5)
  logRr <- c(2)
  seLogRr <- c(2^1023)
  trueLogRr <- c(0.25)
  expect_equal(minLogLikelihoodErrorModelLegacy(theta, logRr, seLogRr, trueLogRr), 99999)
})

# skewNormalLlApproximation
test_that("skewNormalLlApproximation - is infinite sigma", {
  x <- c(1, 2, 3)
  parameters <- list()
  parameters$sigma <- Inf
  expect_equal(skewNormalLlApproximation(x, parameters), c(0, 0, 0))
})
