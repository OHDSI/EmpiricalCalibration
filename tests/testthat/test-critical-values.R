library(testthat)
library(EmpiricalCalibration)

test_that("check computeCvPoisson boundary conditions", {
  expect_error(computeCvPoisson(c(-1, 2)), "should be positive")
  expect_error(computeCvPoisson(c(0, 1, 2), minimumEvents = 1.3), "positive integer")
  expect_error(computeCvPoisson(c(0, 1, 2), minimumEvents = -1), "positive integer")
  expect_error(computeCvPoisson(c(0, 1, 2), alpha = -0.4), "between 0 and 1")
  expect_error(computeCvPoisson(c(0, 1, 2), alpha = 1.3), "between 0 and 1")
  expect_error(computeCvPoisson(c(0, 1, 2), sampleSize = 1.3), "positive integer")
  expect_error(computeCvPoisson(c(0, 1, 2), sampleSize = -1), "positive integer")
})

test_that("check computeCvBinomial boundary conditions", {
  expect_error(computeCvBinomial(c(-1, 2)), "should be positive")
  expect_error(computeCvBinomial(c(0, 1, 2), z = -1.2), "positive")
  expect_error(computeCvBinomial(c(0, 1, 2), z = 1, minimumEvents = 1.3), "positive integer")
  expect_error(computeCvBinomial(c(0, 1, 2), z = 1, minimumEvents = -1), "positive integer")
  expect_error(computeCvBinomial(c(0, 1, 2), z = 1, alpha = -0.4), "between 0 and 1")
  expect_error(computeCvBinomial(c(0, 1, 2), z = 1, alpha = 1.3), "between 0 and 1")
  expect_error(computeCvBinomial(c(0, 1, 2), z = 1, sampleSize = 1.3), "positive integer")
  expect_error(computeCvBinomial(c(0, 1, 2), z = 1, sampleSize = -1), "positive integer")
})

test_that("computeCvPoisson has same output as Sequential", {
  skip_if_not_installed("Sequential")
  # Compare the output of Poisson CV computed using the implementation in this package with output from the implementation in Sequential package

  groupSizes <- rep(1, 10)
  goldStandard <- Sequential::CV.Poisson(
    SampleSize = sum(groupSizes),
    M = 1,
    GroupSizes = groupSizes
  )

  cv <- computeCvPoisson(
    groupSizes = groupSizes,
    minimumEvents = 1
  )
  expect_equal(cv, goldStandard, tolerance = 1e-5, check.attributes = FALSE)
})

test_that("computeCvBinomial has same output as Sequential", {
  skip_if_not_installed("Sequential")
  # Compare the output of Binomial CV computed using the implementation in this package with output from the implementation in Sequential package

  groupSizes <- rep(1, 10)
  z <- 4
  goldStandard <- Sequential::CV.Binomial(
    N = sum(groupSizes),
    z = z,
    M = 1,
    GroupSizes = groupSizes
  )

  cv <- computeCvBinomial(
    groupSizes = groupSizes,
    z = z,
    minimumEvents = 1
  )
  expect_equal(cv, goldStandard$cv, tolerance = 1e-5, check.attributes = FALSE)

  expect_equal(attr(cv, "alpha"), goldStandard$Type_I_Error, tolerance = 1e-3, check.attributes = FALSE)
})
