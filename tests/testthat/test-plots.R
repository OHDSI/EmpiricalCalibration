# We still need to test that plots are saved when that is required

library(testthat)
library(EmpiricalCalibration)

data(sccs)
negatives <- sccs[sccs$groundTruth == 0, ]
positives <- sccs[sccs$groundTruth == 1, ]

null <- fitNull(
  logRr   = negatives$logRr,
  seLogRr = negatives$seLogRr
)

mcmcNull <- fitMcmcNull(
  logRr = negatives$logRr,
  seLogRr = negatives$seLogRr,
  iter = 1e3
)

test_that("output of plotMcmcTrace", {
  expect_is(
    plotMcmcTrace(
      mcmcNull = mcmcNull
    ),
    "ggplot"
  )
})

test_that("output of plotForest", {
  expect_is(
    plotForest(
      logRr = negatives$logRr,
      seLogRr = negatives$seLogRr,
      names = negatives$drugName
    ),
    "ggplot"
  )
})

test_that("output of plotCalibration", {
  expect_warning(
    plotCalibration(
      c(negatives$logRr, Inf),
      c(negatives$seLogRr, 0)
    ),
    "infinite logRr"
  )

  expect_warning(
    plotCalibration(
      c(negatives$logRr, 0),
      c(negatives$seLogRr, Inf)
    ),
    "infinite standard error"
  )

  expect_warning(
    plotCalibration(
      c(negatives$logRr, NA),
      c(negatives$seLogRr, 0)
    ),
    "NA logRr"
  )

  expect_warning(
    plotCalibration(
      c(negatives$logRr, 0),
      c(negatives$seLogRr, NA)
    ),
    "NA standard error"
  )

  expect_is(
    plotCalibration(negatives$logRr, negatives$seLogRr),
    "ggplot"
  )
})

test_that("output of plotCalibrationEffect", {
  expect_error(
    plotCalibrationEffect(
      logRrNegatives   = negatives$logRr,
      seLogRrNegatives = negatives$seLogRr,
      logRrPositives   = positives$logRr,
      seLogRrPositives = positives$seLogRr,
      showCis          = TRUE,
      null             = null
    ),
    "fitMcmcNull"
  )

  expect_is(
    plotCalibrationEffect(
      logRrNegatives   = negatives$logRr,
      seLogRrNegatives = negatives$seLogRr,
      logRrPositives   = positives$logRr,
      seLogRrPositives = positives$seLogRr
    ),
    "ggplot"
  )

  expect_is(
    plotCalibrationEffect(
      logRrNegatives                      = negatives$logRr,
      seLogRrNegatives                    = negatives$seLogRr,
      logRrPositives                      = positives$logRr,
      seLogRrPositives                    = positives$seLogRr,
      showExpectedAbsoluteSystematicError = TRUE
    ),
    "ggplot"
  )

  expect_is(
    plotCalibrationEffect(
      logRrNegatives                      = negatives$logRr,
      seLogRrNegatives                    = negatives$seLogRr,
      logRrPositives                      = positives$logRr,
      seLogRrPositives                    = positives$seLogRr,
      showExpectedAbsoluteSystematicError = TRUE,
      null                                = null
    ),
    "ggplot"
  )

  expect_is(
    plotCalibrationEffect(
      logRrNegatives   = negatives$logRr,
      seLogRrNegatives = negatives$seLogRr,
      logRrPositives   = positives$logRr,
      seLogRrPositives = positives$seLogRr,
      null             = null
    ),
    "ggplot"
  )

  expect_is(
    plotCalibrationEffect(
      logRrNegatives   = negatives$logRr,
      seLogRrNegatives = negatives$seLogRr,
      null             = null
    ),
    "ggplot"
  )

  expect_is(
    plotCalibrationEffect(
      logRrNegatives   = negatives$logRr,
      seLogRrNegatives = negatives$seLogRr,
      logRrPositives   = positives$logRr,
      seLogRrPositives = positives$seLogRr,
      showCis          = TRUE
    ),
    "ggplot"
  )

  # Custom x and y limits
  expect_warning(
    plotCalibrationEffect(
      logRrNegatives = negatives$logRr,
      seLogRrNegatives = negatives$seLogRr,
      logRrPositives = c(-3, -2, 11),
      seLogRrPositives = c(0.1, 1.2, 2.5)
    ),
    regexp = "xLimits"
  )

  expect_warning(
    plotCalibrationEffect(
      logRrNegatives = negatives$logRr,
      seLogRrNegatives = negatives$seLogRr,
      logRrPositives = c(-3, -2, 2),
      seLogRrPositives = c(0.1, 1.2, 2.5)
    ),
    regexp = "yLimits"
  )
})

test_that("output of plotCiCalibration", {
  expect_is(
    plotCiCalibration(
      logRr     = negatives$logRr,
      seLogRr   = negatives$seLogRr,
      trueLogRr = negatives$groundTruth
    ),
    "ggplot"
  )
})

test_that("output of plotCiCalibrationEffect", {
  expect_is(
    plotCiCalibrationEffect(
      logRr = negatives$logRr,
      seLogRr = negatives$seLogRr,
      trueLogRr = negatives$groundTruth
    ),
    "ggplot"
  )
})

test_that("output of plotCiCoverage", {
  expect_is(
    plotCiCoverage(
      logRr = negatives$logRr,
      seLogRr = negatives$seLogRr,
      trueLogRr = negatives$groundTruth
    ),
    "ggplot"
  )
})


test_that("output of plotErrorModel", {
  expect_is(
    plotErrorModel(
      logRr     = negatives$logRr,
      seLogRr   = negatives$seLogRr,
      trueLogRr = negatives$groundTruth
    ),
    "ggplot"
  )
})

test_that("output of plotExpectedType1Error", {
  expect_error(
    plotExpectedType1Error(
      negatives$logRr,
      negatives$seLogRr,
      positives$seLogRr,
      showCis = TRUE,
      null    = null
    ),
    regexp = "fitMcmcNull"
  )

  expect_is(
    plotExpectedType1Error(
      negatives$logRr,
      negatives$seLogRr,
      positives$seLogRr
    ),
    "ggplot"
  )

  expect_is(
    plotExpectedType1Error(
      negatives$logRr,
      negatives$seLogRr,
      positives$seLogRr,
      showCis = TRUE
    ),
    "ggplot"
  )

  expect_is(
    plotExpectedType1Error(
      negatives$logRr,
      negatives$seLogRr,
      positives$seLogRr,
      showEffectSizes = TRUE
    ),
    "gtable"
  )
})
