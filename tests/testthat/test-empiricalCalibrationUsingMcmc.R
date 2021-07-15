library(testthat)
library(EmpiricalCalibration)

data(sccs)
# negatives <- sccs[sccs$groundTruth == 0, ]

test_that("fitMcmcNull requirements", {
  
  # Infinite logRr
  logRr     <- c(0, Inf)
  seLogRr   <- c(1, 0)
  expect_warning(
    fitMcmcNull(
      logRr     = logRr,
      seLogRr   = seLogRr
    ),
    regexp = ".*infinite logRr.*"
  )
  
  # Infinite seLogRr
  logRr     <- c(0, 0)
  seLogRr   <- c(1, Inf)
  expect_warning(
    fitMcmcNull(
      logRr     = logRr,
      seLogRr   = seLogRr
    ),
    regexp = ".*infinite standard error.*"
  )
  
  # NA logRr
  logRr     <- c(0, NA)
  seLogRr   <- c(1, 0)
  expect_warning(
    fitMcmcNull(
      logRr     = logRr,
      seLogRr   = seLogRr
    ),
    regexp = ".*NA logRr.*"
  )
  
  # NA seLogRr
  logRr     <- c(0, 0)
  seLogRr   <- c(1, NA)
  expect_warning(
    fitMcmcNull(
      logRr     = logRr,
      seLogRr   = seLogRr
    ),
    regexp = ".*NA standard error.*"
  )
})
  
