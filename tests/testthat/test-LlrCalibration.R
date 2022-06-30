library(testthat)
library(EmpiricalCalibration)

data(sccs)
negatives <- sccs[sccs$groundTruth == 0, ]
positive <- sccs[sccs$groundTruth == 1, ]

test_that("calibrateLlr normal distribution", {
  null <- fitNull(negatives$logRr, negatives$seLogRr)

  expect_equal(calibrateLlr(null, positive), 0)

  expect_error(calibrateLlr(null, positive, twoSided = TRUE),
    regexp = ".*only one-sided upper LLRs are supported.*"
  )

  expect_error(calibrateLlr(null, positive, upper = FALSE),
    regexp = ".*only one-sided upper LLRs are supported.*"
  )
})

test_that("calibrateLlr grid distribution", {
  set.seed(123)
  mu <- 0.2
  sigma <- 0.2
  data <- simulateControls(n = 50, mean = mu, sd = sigma, trueLogRr = 0, seLogRr = 0.1)
  goldStandardNull <- fitNull(logRr = data$logRr, seLogRr = data$seLogRr)

  point <- seq(log(0.1), log(10), length.out = 1000)
  createGridApproximation <- function(row) {
    return(data.frame(
      point = point,
      value = dnorm(point, mean = row$logRr, sd = row$seLogRr)
    ))
  }
  gridApproximations <- lapply(split(data, 1:nrow(data)), createGridApproximation)
  null <- fitNullNonNormalLl(gridApproximations)
  calibrateLlr(null, gridApproximations)
})

test_that("calibrateLlr non normal distribution", {
  null <- fitNull(negatives$logRr, negatives$seLogRr)
  colnames(positive) <- c("drugName", "mu", "gamma", "sigma")
  expect_equal(calibrateLlr(null, positive), 0.247,
    tolerance = 0.1,
    check.attributes = FALSE
  )


  colnames(positive) <- c("drugName", "mu", "alpha", "sigma")
  expect_equal(calibrateLlr(null, positive), 0.344,
    tolerance = 0.1,
    check.attributes = FALSE
  )
})
