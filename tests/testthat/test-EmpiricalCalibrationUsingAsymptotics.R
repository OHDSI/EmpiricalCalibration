library(testthat)
library(EmpiricalCalibration)

data(sccs)

test_that("calibrateP.null argument requirements", {
  logRr <- c(0.1, 0.3, 0.2)
  seLogRr <- c(0.05, 0.1)
  null <- fitNull(c(1, 2), c(0, 0))
  expect_error(calibrateP(null, logRr, seLogRr),
    regexp = ".*arguments must be of equal length.*"
  )
})

test_that("calibrateOneP argument requirements", {
  logRr <- c(0.1, 0.3, 0.2)
  seLogRr <- c(0.05, 0.1, NA)
  null <- fitNull(c(0, 0), c(0, 0))
  expect_equal(calibrateP(null, logRr, seLogRr), c(0.045500266, 0.002699796, NA))
  expect_equal(calibrateP(null, logRr, seLogRr, twoSided = FALSE), c(0.022750133, 0.001349898, NA))
  expect_equal(round(calibrateP(null, logRr, seLogRr, twoSided = FALSE, upper = FALSE), 2), c(0.98, 1.00, NA))
})


test_that("computeTraditionalP argument options", {
  positive <- sccs[sccs$groundTruth == 1, ]
  expect_equal(computeTraditionalP(positive$logRr, positive$seLogRr), 0)
  expect_equal(computeTraditionalP(positive$logRr, positive$seLogRr, twoSided = FALSE), 0)
  expect_equal(computeTraditionalP(positive$logRr, positive$seLogRr, twoSided = FALSE, upper = FALSE), 1)
})


test_that("fitNullNonNormalLl", {
  negatives <- sccs[sccs$groundTruth == 0, ]
  null <- fitNull(negatives$logRr, negatives$seLogRr)
  expect_equal(fitNullNonNormalLl(negatives), null)
})

test_that("fitNullNonNormalLl using non-normal approximation", {
  skip_on_cran()
  set.seed(123)

  # Test for fitting null using non-normal approximation and asymptotics
  mu <- 0.2
  sigma <- 0.2
  data <- simulateControls(n = 50, mean = mu, sd = sigma, trueLogRr = 0, seLogRr = 0.1)
  goldStandardNull <- fitNull(logRr = data$logRr, seLogRr = data$seLogRr)

  point <- seq(log(0.1), log(10), length.out = 1000)
  createGridApproximation <- function(row) {
    return(data.frame(
      point = point,
      value = dnorm(point, mean = row$logRr, sd = row$seLogRr)
    ))
  }
  gridApproximations <- lapply(split(data, 1:nrow(data)), createGridApproximation)
  null <- fitNullNonNormalLl(gridApproximations)

  expect_equal(null[1],
    goldStandardNull[1],
    tolerance = 0.1,
    scale = 1,
    check.attributes = FALSE
  )
  expect_equal(null[2],
    goldStandardNull[2],
    tolerance = 0.1,
    scale = 1,
    check.attributes = FALSE
  )

  # Test for fitting null using normal approximation and MCMC
  mu <- 0.2
  sigma <- 0.2
  data <- simulateControls(n = 50, mean = mu, sd = sigma, trueLogRr = 0, seLogRr = 0.01)
  null <- fitMcmcNull(logRr = data$logRr, seLogRr = data$seLogRr)
  mcmc <- attr(null, "mcmc")
  lb99Mu <- quantile(mcmc$chain[, 1], 0.005)
  ub99Mu <- quantile(mcmc$chain[, 1], 0.995)
  lb99Sigma <- 1 / sqrt(quantile(mcmc$chain[, 2], 0.995))
  ub99Sigma <- 1 / sqrt(quantile(mcmc$chain[, 2], 0.005))

  expect_lt(lb99Mu, mu)
  expect_gt(ub99Mu, mu)
  expect_lt(lb99Sigma, sigma)
  expect_gt(ub99Sigma, sigma)
})

test_that("fitNullNonNormalLl test for errors and warnings", {
  negatives <- sccs[sccs$groundTruth == 0, ]
  colnames(negatives) <- c("drugName", "mu", "grid", "sigma")
  expect_error(fitNullNonNormalLl(negatives),
    regexp = ".*but not all column names are numeric.*"
  )

  colnames(negatives) <- c("drugName", "mu", "gamma", "sigma")
  negatives$mu[1] <- NA
  expect_warning(fitNullNonNormalLl(negatives),
    regexp = ".*Approximations with NA parameters detected.*"
  )

  colnames(negatives) <- c("drugName", "mu", "alpha", "sigma")
  negatives$mu[1] <- NA
  expect_warning(fitNullNonNormalLl(negatives),
    regexp = ".*Approximations with NA parameters detected.*"
  )

  colnames(negatives) <- c(1, 2, 3, 4)
  null <- fitNullNonNormalLl(negatives)
  expect_equivalent(round(null["mean"], 2), 0)
  expect_equivalent(round(null["sd"], 2), 0.1)
})

test_that("CalibrateP matches computeTraditionalP when mu = sigma = 0", {
  null <- c(
    mean = 0,
    sd = 0
  )
  class(null) <- "null"
  logRr <- .2
  seLogRr <- .2

  expect_equal(
    calibrateP(null, logRr, seLogRr),
    computeTraditionalP(logRr, seLogRr)
  )
})
