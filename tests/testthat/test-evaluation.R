library(testthat)
library(EmpiricalCalibration)

data(sccs)
negatives <- sccs[sccs$groundTruth == 0, ]
negatives <- negatives[1:5, ]

test_that("evaluateCiCalibration requirements", {
  logRr <- c(0, 0)
  seLogRr <- c(1, 0)
  trueLogRr <- c(0, 0)
  strata <- c("strat1", "strat2")
  expect_error(
    evaluateCiCalibration(
      logRr     = logRr,
      seLogRr   = seLogRr,
      trueLogRr = trueLogRr,
      strata    = strata
    ),
    regexp = "factor"
  )

  # Infinite logRr
  logRr <- c(negatives$logRr, Inf)
  seLogRr <- c(negatives$seLogRr, 0)
  trueLogRr <- c(negatives$groundTruth, 0)
  expect_warning(
    evaluateCiCalibration(
      logRr     = logRr,
      seLogRr   = seLogRr,
      trueLogRr = trueLogRr
    ),
    regexp = ".*infinite logRr.*"
  )

  # Infinite seLogRr
  logRr <- c(negatives$logRr, 0)
  seLogRr <- c(negatives$seLogRr, Inf)
  trueLogRr <- c(negatives$groundTruth, 0)
  expect_warning(
    evaluateCiCalibration(
      logRr     = logRr,
      seLogRr   = seLogRr,
      trueLogRr = trueLogRr
    ),
    regexp = ".*infinite standard error.*"
  )

  # logRr is NA
  logRr <- c(negatives$logRr, NA)
  seLogRr <- c(negatives$seLogRr, 0)
  trueLogRr <- c(negatives$groundTruth, 0)
  expect_warning(
    evaluateCiCalibration(
      logRr     = logRr,
      seLogRr   = seLogRr,
      trueLogRr = trueLogRr
    ),
    regexp = ".*NA logRr.*"
  )

  # seLogRr is NA
  logRr <- c(negatives$logRr, 0)
  seLogRr <- c(negatives$seLogRr, NA)
  trueLogRr <- c(negatives$groundTruth, 0)
  expect_warning(
    evaluateCiCalibration(
      logRr     = logRr,
      seLogRr   = seLogRr,
      trueLogRr = trueLogRr
    ),
    regexp = ".*NA standard error.*"
  )
})
