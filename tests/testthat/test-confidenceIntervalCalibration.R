library(testthat)
library(EmpiricalCalibration)

data(sccs)
negatives <- sccs[sccs$groundTruth == 0, ]

test_that("fitSystematicErrorModel requirements", {
  logRr     <- c(0, 0)
  seLogRr   <- c(1, Inf)
  trueLogRr <- c(0, 0)
  
  # Infinite standard error
  expect_warning(
    fitSystematicErrorModel(
      logRr     = logRr,
      seLogRr   = seLogRr,
      trueLogRr = trueLogRr
    ),
    regexp = ".*infinite standard error"
  )
  
  # Infinite logRr
  logRr   <- c(0, Inf)
  seLogRr <- c(1, 0)
  expect_warning(
    fitSystematicErrorModel(
      logRr     = logRr,
      seLogRr   = seLogRr,
      trueLogRr = trueLogRr
    ),
    regexp = ".*infinite logRr"
  )
  
  # seLogRr is NA
  logRr   <- c(0, 0)
  seLogRr <- c(1, NA)
  expect_warning(
    fitSystematicErrorModel(
      logRr     = logRr,
      seLogRr   = seLogRr,
      trueLogRr = trueLogRr
    ),
    regexp = ".*NA standard error.*"
  )
  
  # logRr is NA
  logRr   <- c(0, NA)
  seLogRr <- c(1, 0)
  expect_warning(
    fitSystematicErrorModel(
      logRr     = logRr,
      seLogRr   = seLogRr,
      trueLogRr = trueLogRr
    ),
    regexp = ".*NA logRr.*"
  )
})

test_that("convertNullToErrorModel requirements", {
  null <- 0
  class(null) <- "mcmcNul"
  expect_error(
    convertNullToErrorModel(null = null),
    regexp = ".*type 'null'.*"
  )
})