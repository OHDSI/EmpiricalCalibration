library(testthat)
library(EmpiricalCalibration)

# First three tests check whether functions properly return theoretical values

test_that("computeExpectedAbsoluteSystematicError returns mean value", {
  data(sccs)
  negatives <- sccs[sccs$groundTruth == 0, ]
  null <- fitNull(negatives$logRr, negatives$seLogRr)
  error <- computeExpectedAbsoluteSystematicError(null)
  expect_equal(error, null["mean"], tolerance = 1e-3, check.attributes = FALSE)
})

test_that("computeExpectedAbsoluteSystematicError.null returns zero value", {
  mean <- 0
  sd <- 0
  null <- c(mean, sd)
  error <- computeExpectedAbsoluteSystematicError.null(null)
  expect_equal(error, 0, tolerance = 0, check.attributes = FALSE)
})

test_that("computeExpectedAbsoluteSystematicError.null returns mean value", {
  data(sccs)
  negatives <- sccs[sccs$groundTruth == 0, ]
  null <- fitNull(negatives$logRr, negatives$seLogRr)
  error <- computeExpectedAbsoluteSystematicError.null(null)
  expect_equal(error, null["mean"], tolerance = 1e-3, check.attributes = FALSE)
})

# Below test checks whether function properly returns values
test_that("computeExpectedAbsoluteSystematicError.mcmcNull returns mean", {
  alpha <- 0.05
  data(sccs)
  negatives <- sccs[sccs$groundTruth == 0, ]
  null <- fitMcmcNull(negatives$logRr, negatives$seLogRr)
  chain <- attr(null, "mcmc")$chain
  dist <- apply(chain, 1, function(x) closedFormIntegeralAbsolute(x[1], 1 / sqrt(x[2])))
  result <- quantile(dist, c(0.5, alpha / 2, 1 - (alpha / 2)))

  error <- computeExpectedAbsoluteSystematicError.mcmcNull(null)

  expect_equal(error[1], result[1], tolerance = 0, check.attributes = FALSE)
  expect_equal(error[2], result[2], tolerance = 0, check.attributes = FALSE)
  expect_equal(error[3], result[3], tolerance = 0, check.attributes = FALSE)
})
