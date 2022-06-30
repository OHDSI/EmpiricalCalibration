library(testthat)
library(EmpiricalCalibration)

data(sccs)
negatives <- sccs[sccs$groundTruth == 0, ]

test_that("fitNull argument requirements", {
  logRr <- c(Inf)
  expect_warning(fitNull(logRr, negatives$seLogRr),
    regexp = ".*infinite logRr.*"
  )

  seLogRr <- c(Inf)
  expect_warning(fitNull(negatives$logRr, seLogRr),
    regexp = ".*infinite standard error.*"
  )

  logRr <- c(NA)
  expect_warning(fitNull(logRr, negatives$seLogRr),
    regexp = ".*NA logRr.*"
  )

  seLogRr <- c(NA)
  expect_warning(fitNull(negatives$logRr, seLogRr),
    regexp = ".*NA standard error.*"
  )
})

test_that("fitNull output", {
  null <- fitNull(c(1, 2), c(0, 0))
  expect_equivalent(null["mean"], 1.5, tolerance = 1E-3)

  text <- capture.output(null)
  expect_equal(text[1], "Estimated null distribution")

  # null <- fitNull(negatives$logRr, negatives$seLogRr)
})
