library(testthat)
library(EmpiricalCalibration)

test_that("simulateControls", {
  n <- 100
  mean <- 2
  trueLogRr <- 1
  data <- simulateControls(n = n, mean = mean, sd = 0, seLogRr = 0, trueLogRr = trueLogRr)
  expect_equal(data$logRr, rep(mean + trueLogRr, n))
})

test_that("plotTrueAndObserved", {
  set.seed(123)
  data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueLogRr = log(c(1, 2, 4)))
  fileName <- file.path(tempdir(), "output.png")
  plot <- plotTrueAndObserved(data$logRr, data$seLogRr, data$trueLogRr, fileName = fileName)
  expect_true(file.exists(fileName))
  expect_false(is.null(plot)) # Obviously this could be much better
})
