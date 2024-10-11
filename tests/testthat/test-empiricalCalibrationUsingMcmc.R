library(testthat)
library(EmpiricalCalibration)

data(sccs)
# negatives <- sccs[sccs$groundTruth == 0, ]

test_that("fitMcmcNull requirements", {

  # Infinite logRr
  logRr <- c(0, Inf)
  seLogRr <- c(1, 0)
  expect_warning(
    fitMcmcNull(
      logRr     = logRr,
      seLogRr   = seLogRr
    ),
    regexp = ".*infinite logRr.*"
  )

  # Infinite seLogRr
  logRr <- c(0, 0)
  seLogRr <- c(1, Inf)
  expect_warning(
    fitMcmcNull(
      logRr     = logRr,
      seLogRr   = seLogRr
    ),
    regexp = ".*infinite standard error.*"
  )

  # NA logRr
  logRr <- c(0, NA)
  seLogRr <- c(1, 0)
  expect_warning(
    fitMcmcNull(
      logRr     = logRr,
      seLogRr   = seLogRr
    ),
    regexp = ".*NA logRr.*"
  )

  # NA seLogRr
  logRr <- c(0, 0)
  seLogRr <- c(1, NA)
  expect_warning(
    fitMcmcNull(
      logRr     = logRr,
      seLogRr   = seLogRr
    ),
    regexp = ".*NA standard error.*"
  )
})

test_that("MCMC calibration of p-values returns values close to truth", {
  skip_on_cran()
  set.seed(124)

  negatives <- simulateControls(
    n         = 50,
    mean      = 0.25,
    sd        = 0.16,
    trueLogRr = 0
  )

  # Testing if the fitted MCMC null is correct
  null <- fitMcmcNull(
    logRr   = negatives$logRr,
    seLogRr = negatives$seLogRr,
    iter    = 1e4
  )

  expect_equal(null[1], .25, tolerance = .01)

  # Testing calibrated p-values (lower or two-sided)
  testValue <- calibrateP(
    null     = null,
    logRr    = .2,
    seLogRr  = .2,
    twoSided = FALSE,
    upper    = FALSE
  )

  expect_equal(
    testValue$p,
    pnorm(.2, .25, .16 + .2^2),
    tolerance = .01
  )

  testValue <- calibrateP(
    null     = null,
    logRr    = .2,
    seLogRr  = .2,
    twoSided = TRUE
  )
  target <- min(
    pnorm(.2, .25, .16 + .2^2, lower.tail = TRUE),
    pnorm(.2, .25, .16 + .2^2, lower.tail = FALSE)
  )

  expect_equal(
    testValue$p,
    2 * target,
    tolerance = .01
  )

  # Testing return output of NA logRr or seLogRr
  calibration <- calibrateP(
    null     = null,
    logRr    = NA,
    seLogRr  = .2,
    twoSided = FALSE,
    upper    = FALSE
  )
  testValue <- data.frame(
    p      = calibration$p,
    lb95ci = calibration$lb95ci,
    ub95ci = calibration$ub95ci
  )

  expect_equal(
    testValue,
    data.frame(
      p      = as.numeric(NA),
      lb95ci = as.numeric(NA),
      ub95ci = as.numeric(NA)
    )
  )

  calibration <- calibrateP(
    null     = null,
    logRr    = .2,
    seLogRr  = NA,
    twoSided = FALSE,
    upper    = FALSE
  )
  testValue <- data.frame(
    p      = calibration$p,
    lb95ci = calibration$lb95ci,
    ub95ci = calibration$ub95ci
  )

  expect_equal(
    testValue,
    data.frame(
      p      = as.numeric(NA),
      lb95ci = as.numeric(NA),
      ub95ci = as.numeric(NA)
    )
  )
})

test_that("fitMcmcNull throws no error when all estimates are NA", {
  null <- fitMcmcNull(logRr = c(NA, NA, NA),
                      seLogRr = c(NA, NA, NA))
  print(null)
  
  p <- calibrateP(null, 1, 1)
  expect_true(is.na(p$p))
  
  p <- calibrateP(null, 1, 1, pValueOnly = TRUE)
  expect_true(is.na(p))
  
  model <- convertNullToErrorModel(null)
  expect_true(is.na(model[1]))
  
  ci <- calibrateConfidenceInterval(1, 1, model)
  expect_true(is.na(ci[1]))
  
  ease <- computeExpectedAbsoluteSystematicError(null)
  expect_true(is.na(ease$ease))
})
