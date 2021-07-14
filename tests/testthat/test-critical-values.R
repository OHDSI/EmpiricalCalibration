library(testthat)
library(EmpiricalCalibration)

test_that("computeCvPoisson has same output as Sequential", {
  # Compare the output of the Sequential package with the one implemented in this package
  groupSizes <- rep(1, 10)  
  goldStandard <- Sequential::CV.Poisson(SampleSize = sum(groupSizes),
                                         M = 1,
                                         GroupSizes = groupSizes)
  
  cv <- computeCvPoisson(groupSizes = groupSizes,
                         minimumEvents = 1)
  expect_equal(cv, goldStandard, tolerance = 1e-5, check.attributes = FALSE)
})

test_that("computeCvBinomial has same output as Sequential", {
  groupSizes <- rep(1, 10)  
  z <- 4
  goldStandard <- Sequential::CV.Binomial(N = sum(groupSizes),
                                          z = z,
                                          M = 1,
                                          GroupSizes = groupSizes)
  
  cv <- computeCvBinomial(groupSizes = groupSizes,
                          z = z,
                          minimumEvents = 1)
  expect_equal(cv, goldStandard$cv, tolerance = 1e-5, check.attributes = FALSE)
  
  expect_equal(attr(cv, "alpha"), goldStandard$Type_I_Error, tolerance = 1e-3, check.attributes = FALSE)
})

