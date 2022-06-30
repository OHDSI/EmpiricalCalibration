library(testthat)
library(EmpiricalCalibration)

data(sccs)
negatives <- sccs[sccs$groundTruth == 0, ]

test_that("fitSystematicErrorModel requirements", {
  logRr <- c(0, 0)
  seLogRr <- c(1, Inf)
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
  logRr <- c(0, Inf)
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
  logRr <- c(0, 0)
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
  logRr <- c(0, NA)
  seLogRr <- c(1, 0)
  expect_warning(
    fitSystematicErrorModel(
      logRr     = logRr,
      seLogRr   = seLogRr,
      trueLogRr = trueLogRr
    ),
    regexp = ".*NA logRr.*"
  )

  controls <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueLogRr = log(c(1, 2, 4)))
  model <- fitSystematicErrorModel(controls$logRr, controls$seLogRr, controls$trueLogRr,
    legacy = TRUE, estimateCovarianceMatrix = TRUE
  )
  expect_equal(model[1], 0.254, tolerance = 0.1, check.attributes = FALSE)
})

test_that("convertNullToErrorModel requirements", {
  null <- 0
  class(null) <- "mcmcNul"
  expect_error(
    convertNullToErrorModel(null = null),
    regexp = ".*type 'null'.*"
  )
})

test_that("convertNullToErrorModel", {
  data(sccs)
  negatives <- sccs[sccs$groundTruth == 0, ]
  null <- fitNull(negatives$logRr, negatives$seLogRr)
  model <- convertNullToErrorModel(null)
  positive <- sccs[sccs$groundTruth == 1, ]
  result <- calibrateConfidenceInterval(positive$logRr, positive$seLogRr, model)
  expect_equal(result$logRr, -0.0593, tolerance = 0.1, check.attributes = FALSE)
})
